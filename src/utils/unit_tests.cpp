#include <cassert>
#include <iostream>
using namespace std;

#define GENERIC 

#include "swarm/gpu/pair_calculation.hpp"


template<int nbod>
struct PairTest {
	bool table[nbod][nbod];

	/**
	 * Initialize the table and set the diagonal true because
	 * the diagonal will not be crossed off during the test
	 */
	PairTest() {
		for(int i = 0; i < nbod; i++) 
			for(int j = 0; j < nbod; j++) 
				table[i][j] = (i == j);
	}

	void print_table(){
		cout << " ";
		for(int j = 0; j < nbod; j++) cout << "-";
		cout << endl;
		for(int i = 0; i < nbod; i++) {
			cout << "|";
			for(int j = 0; j < nbod; j++)
				cout << (table[i][j] ? "*" : " ");
			cout << "|\n";
		}
		cout << " ";
		for(int j = 0; j < nbod; j++) cout << "-";
		cout << endl;
	}

	void run(){
		for(int ij = 0; ij < nbod * (nbod-1) /2 ; ij++) {
			int f = first<nbod>(ij), s= second<nbod>(ij);
			assert(f >= 0 && f < nbod); assert(s >= 0 && s < nbod);
			cout << f << ", " << s << endl;

			// Cross of the pair s,f and f,s from the table
			table[s][f] = table[f][s] = true;
		}
	}

	bool pass(){
		bool ret = true;
		for(int i = 0; i < nbod; i++) 
			for(int j = 0; j < nbod; j++)
				ret = ret && table[i][j];
		return ret;
	}

	bool test(){
		run();
		print_table();
		return pass();
	}
};


int main() {
	bool p = 
		PairTest<3>().test() && 
		PairTest<4>().test() && 
		PairTest<5>().test() && 
		PairTest<6>().test() && 
		PairTest<7>().test() && 
		PairTest<8>().test() && 
		PairTest<9>().test() && 
		true;

	cout << "Test " << (p ? "pass" : "failed") << endl;
	return p ? 0 : 1;
}
