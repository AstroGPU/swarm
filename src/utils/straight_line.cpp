#include "swarm/swarm.h"
#include <iostream>
using namespace swarm;
using namespace std;

const int nbod = 3;
const int nsys = 64;

int main(int argc, char* argv[]){
	if(argc <= 1){
		cout << "Usage: straight_line <outputfilename>" << endl;
	}
	const string outputfn = argv[1];
	init(config());
	defaultEnsemble ens = defaultEnsemble::create(nbod,nsys);
	
	for(int s = 0; s < nsys ; s++){
		ens[s].id() = 0;
		ens[s].time() = 0;
		
		// Set all the bodies to zero
		for(int b = 0; b < nbod; b++){
			ens[s][b].mass() = 0;
			for(int c = 0; c < 3; c++)
				ens[s][b][c].pos() = 0, ens[s][b][c].vel() = 0;
		}
		
		// the central one has to pull everything together
		ens[s][0].mass() = 1;
		ens[s][1][0].pos() = -10;
		ens[s][2][0].pos() = +10;

	}
	snapshot::save_text(ens,outputfn);
}
