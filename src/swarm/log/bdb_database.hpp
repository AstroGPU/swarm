#pragma once


#include "log.hpp"
#include <db_cxx.h>


namespace swarm { namespace log {

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
    void openForReading(const std::string& fileName);
    void create(const std::string& fileName);
    void createEmpty(const std::string& fileName);
    
    void put(gpulog::logrecord& lr);

    void close();

private:
    void openInternal(const std::string& fileName, int open_mode);


    DbEnv* env;
    Db primary,
       system_idx,
       time_idx,
       event_idx;

};


} } // close namespace log :: swarm
