#include "../common.hpp"

#include <limits>

#include "bdb_database.hpp"

namespace swarm { namespace log {

using gpulog::logrecord;

const int CACHESIZE = 1024*1024*64 ;


typedef uint32_t sysid_t;
typedef uint8_t evtid_t;
struct pkey_t {
    float time;
    uint32_t system_event_id;

    pkey_t(const float& t,const int& sysid, const int& evid)
        :time(t),system_event_id((uint32_t)evid << 24 | sysid){}

    sysid_t system_id()const{ return (system_event_id & ~ (255 << 24)); }
    evtid_t event_id()const{ return (evtid_t) (system_event_id >> 24); }
};

bool operator <(const pkey_t& a, const pkey_t& b){
    if(a.time == b.time)
        return a.system_event_id < b.system_event_id;
    else
        return a.time < b.time;
}


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
    pkey_t& pkey = *(pkey_t*) key->get_data();
    put_in_dbt(pkey.system_id(), result);
    return 0;
}

int lr_extract_evtid(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
    pkey_t& pkey = *(pkey_t*) key->get_data();
    put_in_dbt(pkey.event_id(), result);
    return 0;
}

int lr_extract_time(Db *secondary, const Dbt *key, const Dbt *data, Dbt *result) {
    pkey_t& pkey = *(pkey_t*) key->get_data();
    put_in_dbt(pkey.time , result);
    return 0;
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
    primary.set_bt_compare(bdb_compare<pkey_t>);
	primary.open(NULL, (baseFileName+".p.db").c_str(), NULL, DB_BTREE, open_mode, 0);

    // Open up the system index database, it has to support
    // duplicates and it is given a smaller cache size
    system_idx.set_cachesize(0,CACHESIZE/4,0);
	system_idx.set_flags(DB_DUP | DB_DUPSORT);
    system_idx.set_bt_compare(bdb_compare<sysid_t>);
    system_idx.set_dup_compare(bdb_compare<pkey_t>);
	system_idx.open(NULL, (baseFileName+".sys.db").c_str(), NULL, DB_BTREE, open_mode , 0);

    // Open up the time index database, it has to support
    // duplicates because our index is not a unique index and
    // it takes a smaller cache size
    time_idx.set_cachesize(0,CACHESIZE/4,0);
	time_idx.set_flags(DB_DUP | DB_DUPSORT);
    time_idx.set_bt_compare(bdb_compare<float>);
    time_idx.set_dup_compare(bdb_compare<pkey_t>);
	time_idx.open(NULL, (baseFileName+".time.db").c_str(), NULL, DB_BTREE, open_mode  , 0);

    event_idx.set_cachesize(0,CACHESIZE/4,0);
	event_idx.set_flags(DB_DUP | DB_DUPSORT);
    event_idx.set_bt_compare(bdb_compare<evtid_t>);
    event_idx.set_dup_compare(bdb_compare<pkey_t>);
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
    
void bdb_database::put(logrecord& lr){
    
    double time; int sys;

	// Based on the implementation in query.cpp
	switch(lr.msgid()){
	case 1: case 2: case 11: case 15: case 16:
		lr >> time >> sys;
        break;

	default:
        sys = 0;
        time = std::numeric_limits<double>::quiet_NaN();
        break;
	}

    pkey_t pkey( (float) time, sys, lr.msgid());
    
    Dbt key(&pkey,sizeof(pkey));
    Dbt data((void*)lr.ptr,lr.len());
    primary.put(NULL,&key,&data,0);
}

void bdb_database::close(){
	primary.close(0);
	system_idx.close(0);
	time_idx.close(0);
	event_idx.close(0);
}

} } // close namespace log :: swarm
