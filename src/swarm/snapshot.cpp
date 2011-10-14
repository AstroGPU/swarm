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

namespace snapshot {

template<class T>
void readfromFILE(FILE*f,T& t,const string& filename = ""){
	if(fread(&t,1,sizeof(t),f) != sizeof(t)){
		throw readfileexception(filename,"File I/O error");
	}
}
template<class T>
void writetoFILE(FILE*f,const T& t,const string& filename = ""){
	if(fwrite(&t,1,sizeof(t),f) != sizeof(t)){
		throw writefileexception(filename,"File I/O error");
	}
}

defaultEnsemble load(const string& filename) throw (readfileexception){
	FILE* f = fopen(filename.c_str(),"rb");

	header h; sys s; body b;

	readfromFILE(f,h,filename);
	hostEnsemble ens = hostEnsemble::create(h.nbod,h.nsys);

	for(int i = 0; i < h.nsys; i++){
		ensemble::SystemRef sr = ens[i];

		readfromFILE(f,s,filename);

		sr.id() = s.id;
		sr.time() = s.time; sr.state() = s.state;

		for(int j = 0; j < h.nbod; j++){
			readfromFILE(f,b,filename);

			sr[j].mass() = b.mass;
			sr[j][0].pos() = b.pos[0];
			sr[j][1].pos() = b.pos[1];
			sr[j][2].pos() = b.pos[2];
			sr[j][0].vel() = b.vel[0];
			sr[j][1].vel() = b.vel[1];
			sr[j][2].vel() = b.vel[2];
			for(int i = 0; i < NUM_PLANET_ATTRIBUTES; i++)
				sr[j].attribute(i) = b.attribute[i];
		}
	}
	fclose(f);
	return ens;
}


const char* DEFAULT_IO_TAG = "SwarmDataFile";
const int CURRENT_IO_VERSION = 1;

defaultEnsemble load_text(const string& filename) throw (readfileexception){
	FILE* f = fopen(filename.c_str(),"r");

	header h; sys s; body b;

	char tag[1024];
	int version = 0;

	fscanf(f, "%s %i\n" , tag, &version);
	if(strcmp( DEFAULT_IO_TAG, tag ) != 0)
		throw readfileexception(filename,"Invalid file, header doesn't match");

	if(version != CURRENT_IO_VERSION )
		throw readfileexception(filename, "Incorrect version");

	fscanf(f,"%i %i\n\n\n",&h.nbod,&h.nsys);
	hostEnsemble ens = hostEnsemble::create(h.nbod,h.nsys);

	for(int i = 0; i < h.nsys; i++){
		ensemble::SystemRef sr = ens[i];

		fscanf(f,"%i %le %i\n", &sr.id(), &sr.time(), &sr.state());

		for(int j = 0; j < h.nbod; j++){
			fscanf(f,"\t%le\n\t%le %le %le\n\t%le %le %le\n\n",
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

void save(defaultEnsemble& ens, const string& filename)  throw (writefileexception){
	FILE* f = fopen(filename.c_str(),"wb");

	header h; sys s; body b;

	h.nsys = ens.nsys(), h.nbod = ens.nbod();

	writetoFILE(f,h,filename);

	for(int i = 0; i < h.nsys; i++){

		ensemble::SystemRef sr = ens[i];
		s.id = sr.id();
		s.time = sr.time(), s.state = sr.state(); 

		writetoFILE(f,s,filename);

		for(int j = 0; j < h.nbod; j++){
			b.pos[0] = sr[j][0].pos();
			b.pos[1] = sr[j][1].pos();
			b.pos[2] = sr[j][2].pos();
			b.vel[0] = sr[j][0].vel();
			b.vel[1] = sr[j][1].vel();
			b.vel[2] = sr[j][2].vel();
			b.mass = sr[j].mass();
			for(int i = 0; i < NUM_PLANET_ATTRIBUTES; i++)
				b.attribute[i] = sr[j].attribute(i);

			writetoFILE(f,b,filename);
		}

	}
	fclose(f);
}

void save_text(defaultEnsemble& ens, const string& filename)  throw (writefileexception){
	FILE* f = fopen(filename.c_str(),"w");


	fprintf(f, "%s %i\n" , DEFAULT_IO_TAG, CURRENT_IO_VERSION );
	fprintf(f,"%i %i\n\n\n", ens.nbod(), ens.nsys() );

	for(int i = 0; i < ens.nsys(); i++){

		ensemble::SystemRef sr = ens[i];
		fprintf(f,"%i %.15le %i\n", sr.id(), sr.time(), sr.state()); 

		for(int j = 0; j < ens.nbod(); j++){
			fprintf(f,"\t%.15le\n\t%.15le %.15le %.15le\n\t%.15le %.15le %.15le\n\n",
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
}
