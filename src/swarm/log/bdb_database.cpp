#include "../common.hpp"

#include <libgen.h>
#include <unistd.h>
#include <limits>
#include <boost/concept_check.hpp>

#include "bdb_database.hpp"

namespace swarm { namespace log {

using gpulog::logrecord;

const int CACHESIZE = 1024*1024*64 ;

const char* fileFormatVersion = "1";
const char* swarmngVersion = "1.1";


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

const int create_mode = 00600;

void bdb_database::openEnv(const std::string& basedir){
    std::cerr << "Initializing BDB environment in `" << basedir << "`" << std::endl;
    env->open(basedir.c_str(),DB_CREATE | DB_INIT_CDB | DB_INIT_MPOOL,create_mode);
}

DbEnv* bdb_database::createDefaultEnv(){
    DbEnv* env = new DbEnv(0);
    env->set_cachesize(0,CACHESIZE,0);
    return env;
}

std::string directory_name(const std::string& s){
    //std::cerr << "dirname for `"<< s << "`";
    char* w = strdup(s.c_str());
    std::string d(dirname(w));
    free(w);
    //std::cerr << " is `" << d << "`" << std::endl;
    return d;
}

std::string base_name(const std::string& s){
    //std::cerr << "basename for `"<< s << "`";
    char* w = strdup(s.c_str());
    std::string d(basename(w));
    free(w);
    //std::cerr << " is `" << d << "`" << std::endl;
    return d;
}

void bdb_database::openInternal(const std::string& pathName, int open_mode){

    openEnv(directory_name(pathName));
    std::string fileName = base_name(pathName);

    const char * fn = fileName.c_str();
    
    metadata.open(NULL, fn, "metadata", DB_BTREE, open_mode, create_mode);

    primary.set_bt_compare(bdb_compare<pkey_t>);
	primary.open(NULL, fn, "primary", DB_BTREE, open_mode, create_mode);

    // Open up the system index database, it has to support
    // duplicates and it is given a smaller cache size
   // system_idx.set_cachesize(0,CACHESIZE,0);
	system_idx.set_flags(DB_DUP | DB_DUPSORT);
    system_idx.set_bt_compare(bdb_compare<sysid_t>);
    system_idx.set_dup_compare(bdb_compare<pkey_t>);
	system_idx.open(NULL, fn, "system_idx", DB_BTREE, open_mode , create_mode);

    // Open up the time index database, it has to support
    // duplicates because our index is not a unique index and
    // it takes a smaller cache size
  //  time_idx.set_cachesize(0,CACHESIZE,0);
	time_idx.set_flags(DB_DUP | DB_DUPSORT);
    time_idx.set_bt_compare(bdb_compare<float>);
    time_idx.set_dup_compare(bdb_compare<pkey_t>);
	time_idx.open(NULL, fn, "time_idx", DB_BTREE, open_mode  , create_mode);

  //  event_idx.set_cachesize(0,CACHESIZE,0);
	event_idx.set_flags(DB_DUP | DB_DUPSORT);
    event_idx.set_bt_compare(bdb_compare<evtid_t>);
    event_idx.set_dup_compare(bdb_compare<pkey_t>);
	event_idx.open(NULL, fn, "event_idx", DB_BTREE, open_mode , create_mode);

    // Associate the primary table with the indices
    // the lr_extract_* is the function that defines
    // the indexing scheme
	primary.associate(NULL, &system_idx,  &lr_extract_sysid, DB_IMMUTABLE_KEY);
	primary.associate(NULL, &time_idx  ,  &lr_extract_time , DB_IMMUTABLE_KEY);
	primary.associate(NULL, &event_idx ,  &lr_extract_evtid, DB_IMMUTABLE_KEY);
}

void bdb_database::openForReading(const std::string& fileName) {
    openInternal(fileName, DB_RDONLY);
    validateVersionInfo();
}

void bdb_database::create(const std::string& fileName){
    openInternal(fileName, DB_CREATE);
    fillVersionInfo();
}

void bdb_database::createEmpty(const std::string& fileName){
    openInternal(fileName, DB_CREATE );
    fillVersionInfo();
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

void bdb_database::addMetaData(const std::string name, const std::string value){
  Dbt key((void *)name.data(),name.size()), data((void*)value.data(), value.size());
  metadata.put(NULL,&key,&data,0);
}
std::string bdb_database::getMetaData(const std::string name) {
  Dbt key((void*)name.data(),name.size()), data;
  data.set_flags(DB_DBT_MALLOC);
  
  metadata.get(NULL,&key,&data,0);
  
  // Create a new string and free the buffer
  size_t n = data.get_size();
  char* ptr = (char*) data.get_data();
  std::string value(data.get_size(),0);
  std::copy(ptr, ptr+n, value.begin());
  free(ptr); data.set_data(0);
  
  return value;
}

void bdb_database::fillVersionInfo() {
  addMetaData("fileFormatVersion", fileFormatVersion);
  addMetaData("swarmngVersion", swarmngVersion);
}
bool bdb_database::validateVersionInfo() {
  return (getMetaData("fileFormatVersion") == fileFormatVersion ) ;
    // && (getMetaData("swarmngVersion") == swarmngVersion );
}

Pprimary_cursor_t bdb_database::primary_cursor(){
    shared_ptr<primary_cursor_t> c(new primary_cursor_t);
    primary.cursor(0,&c->_c,0);
    return c;
}

void bdb_database::close(){
	event_idx.close(0);
	time_idx.close(0);
	system_idx.close(0);
	primary.close(0);
}



void primary_cursor_t::close(){
    _c->close();
}

bool primary_cursor_t::get(pkey_t& key, lrw_t& lr, uint32_t flags){
    Dbt k;
    Dbt d;
    k.set_data(&key);
    d.set_data(lr.ptr);
    k.set_ulen(sizeof(key));
    d.set_ulen(lr.len);
    d.set_flags(DB_DBT_USERMEM);
    k.set_flags(DB_DBT_USERMEM);
    return _c->get(&k,&d,flags) == 0;
}

bool primary_cursor_t::position_at(pkey_t& key,lrw_t& lr){
    Dbt k;
    Dbt d;
    k.set_data(&key);
    d.set_data(lr.ptr);
    k.set_ulen(sizeof(key));
    d.set_ulen(lr.len);
    d.set_flags(DB_DBT_USERMEM);
    k.set_flags(DB_DBT_USERMEM);
    k.set_size(sizeof(key));
    return _c->get(&k,&d,DB_SET_RANGE) == 0;
}


void bdb_database::flush()
{
  primary.sync(0);time_idx.sync(0); event_idx.sync(0); system_idx.sync(0);
}


} } // close namespace log :: swarm
