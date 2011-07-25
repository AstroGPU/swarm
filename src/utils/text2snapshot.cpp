#include "swarm/snapshot.hpp"
#include "utils.hpp"

int main(int argc,char * argv[]){
	string ifn, ofn;
	parse_cmd(argc,argv,ifn,ofn);
	swarm::defaultEnsemble ens = swarm::snapshot::load_text(ifn);
	swarm::snapshot::save(ens,ofn);
}

