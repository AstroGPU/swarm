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

/** \file tutorial_simple.cpp
 *  \brief A simple tutorial for beginners to use the swarm API. 
 *
 */


/*
 *  This is a simple tutorial used in doxygen pages
 *  should go through program2doxygen before it
 *  can be used by doxygen.
 *  
 *
 */
// \page TutorialBeginner Beginner tutorial for using the API
//
// You can write your own scenarios for integrating ensembles and use parts
// and pieces of the Swarm-NG library for your C++ application.
//
// In this tutorial we go through a simple C++ program to show you how you can
// use Swarm-NG library. You can find the source file at src/tutorials/tutorial_simple.cpp
//
//
// We start at the first line of the program. You need to include needed headers to use the library.
//
// Some standard C++ libraries
#include <iostream>
// Include Swarm library headers. You may need to add other headers if you want to use
// More advanced feautures. But for simple ensemble creation and integration the following
// line should be enough
#include "swarm/swarm.h"
// We define a shortcut; because we need to use cudaThreadSynchronize() after
// every GPU call to make sure that we are in-sync with GPU.
#define SYNC cudaThreadSynchronize()

// All the Swarm-NG classes are enclodes in swarm namespace, a C++ tradition. (std is the standard library)
using namespace swarm;
using namespace std;

// We define parameters for integration as constants so we can change them later.
// We set the destination time for 5 revolutions. Each revolution is 2 * pi radians.
const double destination_time = 5 * 2 * M_PI ;
// Swarm uses the configuration data structure to pass most of parameters to creation functions
// Here we put all of our configuration items in a 2D string array.
const char * config_pairs[5][2] = {  
       { "integrator" , "hermite" }
      ,{ "time_step", "0.001" }
      ,{ "log_writer", "null" } 
      ,{ "nsys" , "4000" } 
      ,{ "nbod" , "3" } 
    };

// The easy way to create a config object is from a 2D array containing pairs of strings.
// You can also use \ref swarm::config::load to load a configuration from a file. 
// We use this config object to configure all parts of swarm. Note that not every function
// reads all the configuration items. For example, generate_ensemble only reads "nsys" and "nbod".
config cfg( config_pairs );

// Since our program is short, we can put everything in the main function.
int main(int argc, char* argv[]){

    // First we create our reference ensemble. We use the automatic generation using generate_ensemble function
    // This function creates a very stable ensemble of systems for us.
    defaultEnsemble ref = generate_ensemble(cfg);

    // We have to make a copy of initial conditions. We would like to compare the results to initial conditions
    // after integration
    defaultEnsemble ens = ref.clone();


    // Initialize Swarm library. Basically, it initializes CUDA system and default logging system.
    swarm::init(cfg);

    // Select and create the integrator. While you can create an integrator by calling its constructor directly.
    // It is recommended to use integrator::create since it gives you more flexibility at runtime.
    // The create function looks in the swarm library and finds the integrator you requested from 
    // the list of available plugins. 
    Pintegrator integ = integrator::create(cfg);

    // Now we set-up the integrator for integrating.
    //
    // First set the ensemble. For a GPU integrator, the GPU memory will be allocated.
    integ->set_ensemble(ens);
    // Now set the destination time where we want to stop the integration.
    integ->set_destination_time ( destination_time );
    // Need to synchronize because \ref integrator::set_ensemble may upload data to GPU.
    SYNC;

    // Now that everything is set-up, it is safe to pull the trigger and 
    // call the integrate method on integrator. Note that since we didn't set
    // a logger for the integrator, it will use the default logging system.
    integ->integrate();
    // Need to synchronize because integrate is a GPU call.
    SYNC;

    // Once the integration done, we need to examine the data. The easiest
    // check is to see if the systems have preserved energy.
    //
    // find_max_energy_conservation_error is a utility function that compares
    // two ensemble of systems. We compare our working ensemble ens to the 
    // unchanged ensemble ref.
    double max_deltaE = find_max_energy_conservation_error(ens, ref );
    std::cout << "Max Energy Conservation Error =  " << max_deltaE << std::endl;

    // This concludes the program, at this point you will need to clean up the data and ensembles.
    // But most of Swarm-NG objects are stored in reference-counter pointers and will be automatically
    // deallocated by C++ runtime library.
    return 0;
}
// If you want to do more with Swarm-NG, you may consider looking into
// following classes
//  -  \ref swarm::ensemble  This is an abstarct class for ensembles, there are trivial methods to access the data
//      in the ensemble. You can examine and change the ensemble the way you want.
//  -  \ref swarm::defaultEnsemble  The concrete class which contains memory management for ensembles. Use \ref defaultEnsemble::create
//     to create new ensembles. The arrays will be automatically allocated and de-allocated based on reference counting.
//  -  \ref swarm::snapshot::load_text  Use this utility function to load an ensemble from a text file.
//  -  \ref swarm::snapshot::save_text  Use this utility function to save an ensemble to a text file.
//  -  \ref swarm::log::manager      Most users won't need it because the default log manager is used by all integrators.
//      But if you need multiple log streams, you need to create your own log manager
//      and use it in integrators by \ref integrator::set_log_manager.
//
// To find the complete listing of this tutorial look at \ref src/tutorials/tutorial_simple.cpp in the source code repository.
