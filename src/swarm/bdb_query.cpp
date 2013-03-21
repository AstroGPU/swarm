#include "bdb_query.hpp"
#include "log/bdb_database.hpp"

namespace swarm { namespace query {

using namespace log;

gpulog::logrecord allocateLogrecord(const int len){
    char* data = new char[len];
    gpulog::internal::header* hdr = (gpulog::internal::header*) data;

    hdr->len = len;
    gpulog::logrecord lr(data);
    return lr;
}

void execute_bdb_query(const std::string &dbfile, time_range_t T, sys_range_t sys, body_range_t bod) {
    DbEnv* env(bdb_database::createDefaultEnv());bdb_database db(env);
    db.openForReading(dbfile);

    Pprimary_cursor_t cur = db.primary_cursor();

    pkey_t key(T.first,0,0); 
    lrw_t lrw(20480);

    cur->position_at(key,lrw);
    cur->prev(key,lrw);
    while(cur->next(key,lrw) && key.time <= T.last){
        //std::cerr << key.time << " % " << key.system_id() << " % " << (int)key.event_id() << std::endl;

        gpulog::logrecord lr = lrw.lr();
        output_record(std::cout, lr, bod); std::cout << std::endl;
    }
    std::cout << std::endl;

    cur->close();

    db.close();
    env->close(0);
}

} } // namespace swarm query
