/*************************************************************************
 * Copyright (C) 2011 by Saleh Dindar and the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

#include "snapshot.hpp"

namespace swarm {

template<class T>
void readfromFILE(FILE*f,T& t){
	if(fread(&t,1,sizeof(t),f) != sizeof(t)){
		throw snapshot::readfileexception();
	}
}
template<class T>
void writetoFILE(FILE*f,const T& t){
	if(fwrite(&t,1,sizeof(t),f) != sizeof(t)){
		throw snapshot::writefileexception();
	}
}

defaultEnsemble snapshot::load(const string& filename) throw (readfileexception){
	FILE* f = fopen(filename.c_str(),"rb");

	header h; sys s; body b;

	readfromFILE(f,h);
	hostEnsemble ens = hostEnsemble::create(h.nbod,h.nsys);

	for(int i = 0; i < h.nsys; i++){
		ensemble::SystemRef sr = ens[i];

		readfromFILE(f,s);

		sr.time() = s.time; sr.active() = s.active;

		for(int j = 0; j < h.nbod; j++){
			readfromFILE(f,b);

			sr[j].mass() = b.mass;
			sr[j][0].pos() = b.pos[0];
			sr[j][1].pos() = b.pos[1];
			sr[j][2].pos() = b.pos[2];
			sr[j][0].vel() = b.vel[0];
			sr[j][1].vel() = b.vel[1];
			sr[j][2].vel() = b.vel[2];
		}
	}
	fclose(f);
	return ens;
}


defaultEnsemble snapshot::load_text(const string& filename) throw (readfileexception){
	FILE* f = fopen(filename.c_str(),"r");

	header h; sys s; body b;

	fscanf(f,"%i %i\n\n\n",&h.nbod,&h.nsys);
	hostEnsemble ens = hostEnsemble::create(h.nbod,h.nsys);

	for(int i = 0; i < h.nsys; i++){
		ensemble::SystemRef sr = ens[i];

		fscanf(f,"%lg %li\n", &sr.time(), &sr.active());

		for(int j = 0; j < h.nbod; j++){
			fscanf(f,"\t%lg\n\t%lg %lg %lg\n\t%lg %lg %lg\n\n",
					&sr[j].mass(),
					&sr[j][0].pos(),
					&sr[j][1].pos(),
					&sr[j][2].pos(),
					&sr[j][0].vel(),
					&sr[j][1].vel(),
					&sr[j][2].vel()
				   );
		}
	}
	fclose(f);
	return ens;
}

void snapshot::save(defaultEnsemble& ens, const string& filename)  throw (writefileexception){
	FILE* f = fopen(filename.c_str(),"wb");

	header h; sys s; body b;

	h.nsys = ens.nsys(), h.nbod = ens.nbod();

	writetoFILE(f,h);

	for(int i = 0; i < h.nsys; i++){

		ensemble::SystemRef sr = ens[i];
		s.time = sr.time(), s.active = sr.active(); 

		writetoFILE(f,s);

		for(int j = 0; j < h.nbod; j++){
			b.pos[0] = sr[j][0].pos();
			b.pos[1] = sr[j][1].pos();
			b.pos[2] = sr[j][2].pos();
			b.vel[0] = sr[j][0].vel();
			b.vel[1] = sr[j][1].vel();
			b.vel[2] = sr[j][2].vel();
			b.mass = sr[j].mass();

			writetoFILE(f,b);
		}

	}
	fclose(f);
}

void snapshot::save_text(defaultEnsemble& ens, const string& filename)  throw (writefileexception){
	FILE* f = fopen(filename.c_str(),"w");

	fprintf(f,"%i %i\n\n\n", ens.nbod(), ens.nsys() );

	for(int i = 0; i < ens.nsys(); i++){

		ensemble::SystemRef sr = ens[i];
		fprintf(f,"%lg %li\n", sr.time(), sr.active()); 

		for(int j = 0; j < ens.nbod(); j++){
			fprintf(f,"\t%lg\n\t%lg %lg %lg\n\t%lg %lg %lg\n\n",
					sr[j].mass(),
					sr[j][0].pos(),
					sr[j][1].pos(),
					sr[j][2].pos(),
					sr[j][0].vel(),
					sr[j][1].vel(),
					sr[j][2].vel()
				   );
		}

		fprintf(f,"\n");

	}
	fclose(f);
}


}
