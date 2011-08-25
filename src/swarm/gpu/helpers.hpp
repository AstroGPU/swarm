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

#pragma once


template<int Begin, int End, int Step = 1>
struct Unroller {
	template<typename Action>
	__device__	static void step(const Action& action) {
			action(Begin);
			Unroller<Begin+Step, End, Step>::step(action);
		}
};

template<int End, int Step>
struct Unroller<End, End, Step> {
	template<typename Action>
	__device__	static void step(const Action& action) {
		}
};


template<template<int> class T,int N,int MAXN,class B,class P>
struct choose {
B operator ()(const int& n, const P& x){
	if(n == N)
		return T<N>::choose(x);
	else if(n <= MAXN && n > N)
		return choose<T,(N<MAXN ? N+1 : MAXN),MAXN,B,P>()(n,x);
	else
		return B();
}
};

namespace swarm {


template<int i>
struct params_t {
	const static int n = i;
};

template<class implementation,class T>
__global__ void generic_kernel(implementation* integ,T a) {
	integ->kernel(a);
}


template< class implementation, class T>
void launch_template(implementation* integ, implementation* gpu_integ, T a)
{
	if(integ->get_ensemble().nbod() == T::n) 
		generic_kernel<<<integ->gridDim(), integ->threadDim(), integ->shmemSize() >>>(gpu_integ,a);

}

template<int N>
struct launch_template_choose {
	template<class integ_pair>
	static void choose(integ_pair p){
		params_t<N> a;
		generic_kernel<<<p.first->gridDim(), p.first->threadDim(), p.first->shmemSize() >>>(p.second,a);
	}
};

template<class implementation>
void launch_templatized_integrator(implementation* integ){

	if(integ->get_ensemble().nbod() <= MAX_NBODIES){
		implementation* gpu_integ;
		cudaErrCheck ( cudaMalloc(&gpu_integ,sizeof(implementation)) );
		cudaErrCheck ( cudaMemcpy(gpu_integ,integ,sizeof(implementation),cudaMemcpyHostToDevice) );

		typedef std::pair<implementation*,implementation*> integ_pair ;
		integ_pair p ( integ, gpu_integ );
		int nbod = integ->get_ensemble().nbod();

		choose< launch_template_choose, 3, MAX_NBODIES, void, integ_pair > c;
			c( nbod, p );

		cudaFree(integ);
	} else {
		char b[100];
		snprintf(b,100,"Invalid number of bodies. (only up to %d bodies per system)",MAX_NBODIES);
		ERROR(b);
	}

}


	
}
