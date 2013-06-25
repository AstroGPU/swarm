#pragma once


#include "log.hpp"
#include <db_cxx.h>


namespace swarm { namespace log {

struct primary_cursor_t;

class bdb_database {

public:
    bdb_database(DbEnv* e):
        env(e),
        primary(e, 0),
        system_idx(e, 0),
        time_idx(e, 0),
        event_idx(e, 0)
    {}

    static DbEnv* createDefaultEnv();
    void openEnv(const std::string& basedir);
    void openForReading(const std::string& fileName);
    void create(const std::string& fileName);
    void createEmpty(const std::string& fileName);
    
    void put(gpulog::logrecord& lr);

    void close(); 

    shared_ptr<primary_cursor_t> primary_cursor();

private:
    void openInternal(const std::string& fileName, int open_mode);


    DbEnv* env;
    Db primary,
       system_idx,
       time_idx,
       event_idx;

};

typedef uint32_t sysid_t;
typedef uint8_t evtid_t;
struct pkey_t {
    float time;
    uint32_t system_event_id;

    pkey_t(const float& t = 0.0,const int& sysid = 0, const int& evid = 0)
        :time(t),system_event_id((uint32_t)evid << 24 | sysid){}

    sysid_t system_id()const{ return (system_event_id & ~ (255 << 24)); }
    evtid_t event_id()const{ return (evtid_t) (system_event_id >> 24); }
};

struct lrw_t {
    char* ptr;
    int len;
    gpulog::logrecord lr() { return gpulog::logrecord(ptr); }
    lrw_t(const int len):ptr(new char[len]),len(len){
        gpulog::internal::header* hdr = (gpulog::internal::header*) ptr;
        hdr->len = 3; hdr->msgid = -1;
    }
};

struct primary_cursor_t {
    Dbc* _c;
    void close();
  /*!
   * The logrecord should be a valid buffer, the function currently does not
   * re-allocate the buffer
   */
    bool get(pkey_t& key,lrw_t& lr, uint32_t flags);
    bool first(pkey_t& key,lrw_t& lr){ return get(key,lr,DB_FIRST); }
    bool position_at(pkey_t& key,lrw_t& lr);
    bool prev(pkey_t& key,lrw_t& lr){ return get(key,lr,DB_PREV); }
    bool next(pkey_t& key,lrw_t& lr){ return get(key,lr,DB_NEXT); }
    bool current(pkey_t& key,lrw_t& lr){ return get(key,lr,DB_CURRENT); }
};
typedef shared_ptr<primary_cursor_t> Pprimary_cursor_t;


} } // close namespace log :: swarm
