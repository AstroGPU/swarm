#include "../common.hpp"

#include <limits>

#include "bdb_database.hpp"

namespace swarm { namespace log {

using gpulog::logrecord;

const int CACHESIZE = 1024*1024*64 ;


struct logdb_primary_key {
    double time;
    int system_id;
    int event_id;
    idx_t recno;
};


/**
 * Helper function to put any constant size
 * type into a Dbt struct. the flag DB_DBT_APPMALLOC
 * hints berkeley db that the data is allocated by
 * the application
 */
template<typename T>
void put_in_dbt(const T& t, Dbt* data){
	data->set_flags(DB_DBT_APPMALLOC);
	data->set_size(sizeof(T));
	data->set_data(new T(t));
}


bool extract_from_ptr(void* ptr, size_t size, double& time, int& sys){
	logrecord l((char*)ptr);
	assert(l.len() == size);

	// Based on the implementation in query.cpp
	switch(l.msgid()){
	case 1: case 2: case 11: case 15: case 16:
		l >> time >> sys;
		return true;

	default:
        sys = -1;
        time = std::numeric_limits<double>::quiet_NaN();
		return false;
	}

}

/**
 * Extract the system ID from a log record
 *
 * @param secondary: The secondary parameter is the database handle for the secondary.
 * @param key      : The key parameter is a Dbt referencing the primary key.
 * @param data     : The data parameter is a Dbt referencing the primary data item.
 * @param result   : The result parameter is a zeroed Dbt in which the callback function
 *                   should fill in data and size fields that describe the secondary key or keys.
 */
int lr_extract_sysid(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
    logdb_primary_key& pkey = *(logdb_primary_key*) key->get_data();
    put_in_dbt(pkey.system_id, result);
    return 0;
}

int lr_extract_evtid(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
    logdb_primary_key& pkey = *(logdb_primary_key*) key->get_data();
    put_in_dbt(pkey.event_id, result);
    return 0;
}

int lr_extract_time(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
    logdb_primary_key& pkey = *(logdb_primary_key*) key->get_data();
    put_in_dbt(pkey.time , result);
    return 0;
}

/**
 * The primary key is a struct
 * we just convert the key into the struct
 * and perform lexicographical compare in order
 * of time , system_id and event_id
 */
int compare_logdb_primary_key(DB* db, const DBT *k1, const DBT* k2){
    const size_t len = sizeof(logdb_primary_key);
    if(k1->size < k2->size)
        return -1;
    else if(k1->size > k2->size)
        return 1;
    else{
        if( (k1->size == len) && (k2->size == len) ) {

            logdb_primary_key& a = *(logdb_primary_key*)(k1->data);
            logdb_primary_key& b = *(logdb_primary_key*)(k2->data);

            if(a.time < b.time) return -1;
            else if(a.time > b.time) return 1;
            else if(a.system_id < b.system_id) return -1;
            else if(a.system_id > b.system_id) return 1;
            else if(a.event_id < b.event_id) return -1;
            else if(a.event_id > b.event_id) return 1;
            else if(a.recno < b.recno) return -1;
            else if(a.recno > b.recno) return 1;
            else return 0;
        }else{
            return 0;
        }
    }
}


/**
 * Make a comparison function for a C++ type T that supports "<" operator
 */
template<typename T>
int bdb_compare(DB* db, const DBT *k1, const DBT* k2){
    if(k1->size < k2->size)
        return -1;
    else if(k1->size > k2->size)
        return 1;
    else{
        if( (k1->size == sizeof(T)) && (k2->size == sizeof(T)) ) {
            T& a = *(T*)(k1->data);
            T& b = *(T*)(k2->data);
            if(a < b) return -1;
            else if(b < a) return 1;
            else return 0;
        }else{
            return 0;
        }
    }
}



void bdb_database::openInternal(const std::string& baseFileName, int open_mode){

    primary.set_cachesize(0,CACHESIZE,0);
    primary.set_bt_compare(compare_logdb_primary_key);
	primary.open(NULL, (baseFileName+".p.db").c_str(), NULL, DB_BTREE, open_mode, 0);

    // Open up the system index database, it has to support
    // duplicates and it is given a smaller cache size
    system_idx.set_cachesize(0,CACHESIZE/4,0);
	system_idx.set_flags(DB_DUP | DB_DUPSORT);
    system_idx.set_bt_compare(bdb_compare<int>);
    system_idx.set_dup_compare(compare_logdb_primary_key);
	system_idx.open(NULL, (baseFileName+".sys.db").c_str(), NULL, DB_BTREE, open_mode , 0);

    // Open up the time index database, it has to support
    // duplicates because our index is not a unique index and
    // it takes a smaller cache size
    time_idx.set_cachesize(0,CACHESIZE/4,0);
	time_idx.set_flags(DB_DUP | DB_DUPSORT);
    time_idx.set_bt_compare(bdb_compare<double>);
    time_idx.set_dup_compare(compare_logdb_primary_key);
	time_idx.open(NULL, (baseFileName+".time.db").c_str(), NULL, DB_BTREE, open_mode  , 0);

    event_idx.set_cachesize(0,CACHESIZE/4,0);
	event_idx.set_flags(DB_DUP | DB_DUPSORT);
    event_idx.set_bt_compare(bdb_compare<int>);
    event_idx.set_dup_compare(compare_logdb_primary_key);
	event_idx.open(NULL, (baseFileName+".evt.db").c_str(), NULL, DB_BTREE, open_mode , 0);

    // Associate the primary table with the indices
    // the lr_extract_* is the function that defines
    // the indexing scheme
	primary.associate(NULL, &system_idx,  &lr_extract_sysid, DB_IMMUTABLE_KEY);
	primary.associate(NULL, &time_idx  ,  &lr_extract_time , DB_IMMUTABLE_KEY);
	primary.associate(NULL, &event_idx ,  &lr_extract_evtid, DB_IMMUTABLE_KEY);
}

void bdb_database::openForReading(const std::string& baseFileName) {
    openInternal(baseFileName, DB_RDONLY);
}

void bdb_database::create(const std::string& baseFileName){
    openInternal(baseFileName, DB_CREATE);
}

void bdb_database::createEmpty(const std::string& baseFileName){
    openInternal(baseFileName, DB_CREATE | DB_TRUNCATE );
}
    
void bdb_database::put(const logrecord& lr, const idx_t& recno){
    logdb_primary_key pkey;
    Dbt key(&pkey,sizeof(pkey)),data;

    // form the primary key
    pkey.event_id = lr.msgid();
    extract_from_ptr((void*)lr.ptr, lr.len(), pkey.time, pkey.system_id);
    pkey.recno = recno;

    data.set_data((void*)lr.ptr);
    data.set_size(lr.len());
    primary.put(NULL,&key,&data,0);
}

void bdb_database::close(){
	primary.close(0);
	system_idx.close(0);
	time_idx.close(0);
	event_idx.close(0);
}

} } // close namespace log :: swarm
