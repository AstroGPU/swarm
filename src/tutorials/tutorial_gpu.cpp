/*************************************************************************
 * Copyright (C) 2009-2010 by Eric Ford & the Swarm-NG Development Team  *
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

/** \file tutorial_gpu.cpp
 *  \brief A tutorial for using GPU integrators. 
 *
 */


/* This is a tutorial used in doxygen pages
 * should go through program2doxgyen before
 * it can be used by doxygen.
 *
 */
// \page TutorialGPU Advanced tutorial for using GPU integrators
//
//  The basic API gives you control for selecting integrators and
//  using them. If you need to have more control over the integration
//  process and how GPU memory is managed, you need to use more
//  advanced API. This tutorial gives demonstrates how to use
//  API for managing GPU integrators and ensembles.
//
//  Basic set up: include the headers needed and others.
//
#include <iostream>
#include "swarm/swarm.h"
#include "swarm/snapshot.hpp"
#define SYNC cudaThreadSynchronize()
using namespace swarm;
using namespace std;

/// Main program
int main(int argc, char* argv[]){


	// This is the easy way to load config from a file our config file
	// will include parameters to set up a GPU integrator
	//
	config cfg = config::load( "myconfig.cfg" );

	// We load initial conditions from an textual input file
	defaultEnsemble ref = snapshot::load_text( "input.txt" );

	// Make a working copy of the initial conditions
    defaultEnsemble ens = ref.clone();

    // Initialize Swarm library. Basically, it initializes CUDA system and default logging system.
    swarm::init(cfg);

	// Select and initialize the integrator using the config. We assume
	// that the integrator specified is a GPU integrator. Note that we 
	// cannot easily check if the integrator is a GPU integrator. The
	// application will fail if the integrator specified is not a GPU
	// integrator.
	gpu::Pintegrator integ = boost::dynamic_pointer_cast<gpu::integrator>(integrator::create(cfg));

	// Since we are managing GPU memory, we need to make a copy of the ensemble on the GPU. For GPU allocated ensembles we 
	// need to use an object of type deviceEnsemble. Only
	// basic properties of a deviceEnsemble is accessible on
	// the host like nbod() and nsys(). Note that you cannot
	// read the data in a deviceEnsemble on the host.
	// The following line uses cudaMalloc to allocate memory
	// for an ensemble of the same size as ens on the GPU.
	// The uses some cudaMemcpy calls to copy the ensemble
	// to the GPU.
	deviceEnsemble device_ens = ens.cloneTo<deviceEnsemble>();

	// To set the ensemble for the integrator we need to use
	// a different version of integrator::set_ensemble that
	// takes two parameters. The simple version of set_ensemble
	// makes a copy of the given ensemble on the GPU. 
	// 
	// We need to provide the host and device ensemble. 
	// The integrator will keep a reference to the 
	// ensembles.
	integ->set_ensemble( ens, device_ens );

	// This time we do more clever trick for setting the 
	// destination time. If we set the destination time
	// to a predefined value and the ensemble is already
	// past that time, nothing is going to happen. Instead
	// we take the median of times of the systems and
	// add a predefined value to that.
	//
	const double integration_time = 5 * 2 * M_PI ;
	const double begin_time = ens.time_ranges().median ;
	const double destination_time = begin_time + integration_time;
	integ->set_destination_time ( destination_time );

	SYNC;

	// Now we can launch the integration many times 
	// without downloading the results. The 
	// systems won't go past the destination time
	// but we have to integrate many times.
	for(int i= 0; i < 10; i++ ){
		integ->core_integrate();
	}

	// Before we can use the results we have to 
	// download the results. The following call
	// translates to some cudaMemcpy calls
	device_ens.copyTo( ens );

	// TODO: Do something better with the data 
    // find_max_energy_conservation_error is a utility function that compares
    // two ensemble of systems. We compare our working ensemble ens to the 
    // unchanged ensemble ref.
    double max_deltaE = find_max_energy_conservation_error(ens, ref );
    std::cout << "Max Energy Conservation Error =  " << max_deltaE << std::endl;

}
// To find the complete listing of this tutorial look at \ref src/tutorials/tutorial_gpu.cpp in the source code repository.
