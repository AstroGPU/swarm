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



/**
 * \brief Template helper for unrolling of loops
 * This template helper is used for unrolling loops
 * This and the template specialization below provide
 * a pattern matching recursive function that is evaluated
 * at compile time and generates code for unrolling a
 * function
 *
 * \param  Begin  the index to start from
 * \param  End    the end index (non-inclusive)
 * \param  Step   the step parameter (defaults to 1)
 *
 * Usage:
\code
Unroller< 2, 8, 2 >::step( action );
\endcode
 * This generates a code like:
\code
action( 2 );
action( 4 );
action( 6 );
\endcode;
 *
 */
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
	__device__	static void step(const Action& action) { }
};


/** \brief Template object function to choose the appropriate
 * instantiation of a function at compile time
 * 
 * Suppose that you define T as:
\code
template<int N> struct T {
	B choose(P p) { 
	  // do something with p
	  return B( ... );
	}
}
\endcode
 *  
 *  and you want to instantiate T for a finite set of values of 
 *  N like 1,2,3,4. and you want to call the appropriate T
 *  at runtime based on n that is provided at runtime.
 *  You will write:
 *
 *  B b = choose< T, 1, 4, B, P>( n, p );
 *
 *  This will try to instantiate T for values 1 to 4 and calls
 *  T::choose with p for appropriate value of n at runtime.
 *  Note that if n is not in range then it returns B().
 *
 */
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


//! This is a wrapper for a compile time integere value.
//! Because CUDA chokes when encountered with integer template
//! values
template<int i>
struct compile_time_params_t {
	const static int n = i;
};


//! Implementation of the generic_kernel. This is used so we can have
//! arbitrary kernels that are template based and member functions of
//! some class.
template<class implementation,class T>
__global__ void generic_kernel(implementation* integ,T compile_time_param) {
	integ->kernel(compile_time_param);
}

/**
 * 
 */
template< class implementation, class T>
void launch_template(implementation* integ, implementation* gpu_integ, T compile_time_param)
{
	if(integ->get_ensemble().nbod() == T::n) 
		generic_kernel<<<integ->gridDim(), integ->threadDim(), integ->shmemSize() >>>(gpu_integ,compile_time_param);
	else
	  ERROR("Error launching kernel.  Active ensemble has " + inttostr(integ->get_ensemble().nbod()) + " bodies per system.\n");
		
}


/**
 * \brief structure crafted to be used with choose template. 
 * 	 This where the actual kernel is launched. Although all the
 * 	 parameters come from the integ_pair value that is passed in.
 * 	 For more information c.f. choose and launch_templatized_integrator.
 */
template<int N>
struct launch_template_choose {
	template<class integ_pair>
	static void choose(integ_pair p){
		compile_time_params_t<N> compile_time_param;
//		$PRINT( "Launch Kernel " << p.first->gridDim().x << "x" << p.first->gridDim().y 
//				<< "=" << (p.first->gridDim().x * p.first->gridDim().y) <<  ",  "
//				<< p.first->threadDim().x << "x" << p.first->threadDim().y 
//				<< "=" << (p.first->threadDim().x * p.first->threadDim().y) <<  ",  "
//				<< ", " << p.first->shmemSize() );
		generic_kernel<<<p.first->gridDim(), p.first->threadDim(), p.first->shmemSize() >>>(p.second,compile_time_param);
	}
};


/** \brief Global interface for launching a templatized integrator.
 *
 * The passed pointer is expected to be a gpu::integrator. The 
 * function inspects the ensemble of the integrator and launches
 * the templatized version of the integrator that matches the
 * number of bodies.
 *
 * The templates are intantiated for number of bodies ranging
 * from 3 to MAX_NBODIES (defined in config.h.in and adjustable using
 * CMake configurations). 
 * If the number of bodies is not in range, then an error is raised.
 * 
 */
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
