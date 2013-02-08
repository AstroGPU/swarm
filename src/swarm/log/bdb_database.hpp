#pragma once


#include "log.hpp"
#include <db_cxx.h>


namespace swarm { namespace log {

typedef long long idx_t; 

class bdb_database {

public:
    bdb_database():
        primary(NULL, 0),
        system_idx(NULL, 0),
        time_idx(NULL, 0),
        event_idx(NULL, 0)
    {}

    void openForReading(const std::string& baseFileName);
    void create(const std::string& baseFileName);
    void createEmpty(const std::string& baseFileName);
    
    void put(const gpulog::logrecord& lr, const idx_t& recno = 0);

    void close();

private:
    void openInternal(const std::string& baseFileName, int open_mode);


    Db primary,
       system_idx,
       time_idx,
       event_idx;

};


} } // close namespace log :: swarm
