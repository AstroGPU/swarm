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
/*! \file integrator.hpp
 *   \brief Definition of class interface for generic and GPU integrator class
 *
 *   This is the only file that external applications need to include
 *   if they need to use integrators
 *
 */

#pragma once

#include "types/ensemble.hpp"
#include "types/config.hpp"
#include "log/logmanager.hpp"


namespace swarm {

/*! Interface class for all integrators.
 *   This class sets out the general functions to 
 *   use an integrator for external applications.
 *
 *   All CPU integrators should be based on this class. It contains
 *   all the variables required for basic integrations. User variables
 *   can be added in derived instances.
 *
 *   GPU integrators should be based on gpu::integrator.
 *
 *   To add your integrator to swarm library look at plugin development guide.
 *
 */
class integrator {	
	public:
	
	//! Default value for maximum number of iterations. c.f. \ref _max_iterations
	const static int _default_max_iterations;
	//! Default value for maximum number of attempts. c.f. \ref _max_attempts
	const static int _default_max_attempts;
	protected:

	//! ensemble to be integrated
	defaultEnsemble _ens;          

	//! the timestamp at which the integration should finish
	double _destination_time;

	//! log manager for log output, set to manager::default() by default
	log::Pmanager _logman;

	//! log object, obtained from log manager
        //  WARNING: Why do we store a raw pointer here?  
        // Any reason not to just use host_log stored in _logman?
	gpulog::host_log* _log;

	//! Maximum number of iterations per each integration step. c.f. \ref integrate for usage
	int _max_iterations;
	//! Maximum number of attempts to complete the integration. c.f. \ref integrate for usage
	int _max_attempts;

	//! Integrater implementation provided by derived instance
	virtual void launch_integrator() = 0 ;

	public:
	//! Inetgrator class should be configurable. 
	//! Derived instances should also have a constructor with similar signature 
	//! and pass on the config parameter.
	integrator(const config &cfg);

	/*! Interfaces function to integrate, for use by general user.
	 *  Calls launch_integrator() several times up to
	 *  the value specified in \ref _max_attempts until all the
	 *  systems are inactive. Each call can go through
	 *  at most a limited number of iterations specified in
	 *  \ref _max_iterations .
	 *
	 *  To set the parameters for integration use set_ensemble(ens),
	 *  set_destination_time(t), set_log_manager(l) 
	 */
	virtual void integrate();

        //! Flush the host log
        virtual void flush_log() {
	  _logman->flush();	  
	}

	//! Access the ensemble subject to integration
	virtual defaultEnsemble& get_ensemble() {
		return _ens;
	}

	//! Set the ensemble subject to integration
	virtual void set_ensemble(defaultEnsemble& ens) {
		_ens = ens;
	}

	//! Set the time marker to end the integration
	virtual void set_destination_time(const double& destination_time) {
		_destination_time = destination_time;
	}

	/*! Loads an integrator using the plugin system. 
	 * value of cfg["integrator"] is used to identify the 
	 * plugin to be instantiated. The integrator plugin
	 * usually requires additional configuration.
	 *
	 * The "integrator_" suffix is automatically added to
	 * the the value of cfg["integrator"].
	 * Example: To load hermite integrator set cfg["integrator"]="hermite"
	 * and then the plugin "integrator_hermite" will be 
	 * loaded and instantiated.
	 *
	 * If cfg["integrator"] contains an invalid value then
	 * plugin_not_found exception is raised.
	 *
	 * to see the list of plugins available use swarm -P.
	 */
	static shared_ptr<integrator> create(const config &cfg);

	//! accessor function to set the manager for log output
	virtual void set_log_manager(log::Pmanager& l);
        virtual gpulog::host_log* get_host_log();

	virtual void set_max_iterations( const int& mi ) { _max_iterations = mi; }
	virtual void set_max_attempts( const int& ma ) { _max_attempts = ma; }

};
typedef shared_ptr<integrator> Pintegrator;


//! Helper function to calculate number of systems with SYSTEM_ACTIVE flag
int number_of_active_systems(defaultEnsemble ens) ;


/*! GPU-based integrators and other GPU tools
 *
 *   All GPU integrators are containted within this namespace.
 *   Other helper structures and functions like Gravitation 
 *   are also in this namespace.
 *
 */
namespace gpu {

/*! Interface class for all GPU based integrators.
 *   It is very similar to its base class integrator. But it is adapted
 *   for GPU integrators.
 *
 *   This class takes care of GPU logs, uploading and downloding ensemble 
 *   from GPU. Derived instances need to override the pure virtual method
 *   \ref launch_integrator with the CUDA kernel launch code.
 *
 *   Because CUDA is not capable of dereferencing virtual methods, this class
 *   cannot launch kernels. However, template classes for launching kernels
 *   are included in \ref gpu/helpers.hpp.
 */
class integrator : public swarm::integrator {
	typedef swarm::integrator Base;
	protected:
	//! Copy of ensemble on host (kept in sync with \ref _dens).
	hostEnsemble& _hens;
	//! Copy of ensemble on [GPU] device (kept in sync with \ref _hens).
	//! \TODO: use cux auto ptr to make sure we don't have memory leaks
	deviceEnsemble _dens;

	//! GPU log object obtained from log manager.
	gpulog::device_log* _log;

	public: 

	//! Pass on constructor
	integrator(const config &cfg);

	/*! Interfaces function to integrate, for use by general user.
	 *  Launches the GPU integrator kernel several times up to
	 *  the value specified in \ref _max_attempts until all the
	 *  systems are inactive. Each kernel call can go through
	 *  at most a limited number of iterations specified in
	 *  \ref _max_iterations .
	 *
	 *  To set the parameters for integration use set_ensemble(ens),
	 *  set_destination_time(t), set_log_manager(l) 
	 */
	virtual void integrate();

	/** To integrate without any bookkeeping
	 *
	 * This function just executes the integration implementation
	 * The callee has to take care of systems running and
	 * stop conditions and memory management
	 *
	 */
	virtual void core_integrate() {
		launch_integrator();
		flush_log();
	}

	//! Read the GPU log object from log manager and set it
	virtual void set_log_manager(log::Pmanager& l);

	//! Set the GPU log object used for loggin output
	void set_log(gpulog::device_log* log) { _log = log; }

	//! Get the GPU log object used for logging output
        gpulog::device_log* get_device_log() { return _log; }

	//! Set the ensemble, only provide an ensemble on host.
	//! This cals automatically creates a copy of ensemble 
	//! on GPU and keeps it in sync.
	//! Refer to set_ensemble(host_ens,device_ens) if you wish
	//! to manage the GPU memory allocation.
	void set_ensemble(defaultEnsemble& ens) {
		_hens = ens;
		_dens = _hens.cloneTo<deviceEnsemble>();
	}

	/*! Set two host and device ensembles, two ensembles should match.
	 *  Be careful: all data on device_ens will be overwritten when you
	 *  call integrate(). 
	 *  Only use this function if you want to save memory when integrating
	 *  the same ensemble on different integrators.
	 *  Example:
	 *  \code
	 *  defaultEnsemble hens = snapshot::load(...);
	 *  deviceEnsemble  dens = hens.cloneTo<deviceEnsemble>();
	 *
	 *  i1->set_ensemble(hens,dens);
	 *  i2->set_ensemble(hens,dens);
	 *
	 *  i1->integrate();
	 *  i2->integrate();
	 *  i1->integrate();
	 *   .
	 *   .
	 *   .
	 *  \endcode
	 *
	 */
	void set_ensemble(defaultEnsemble& host_ens, deviceEnsemble& device_ens) {
		_hens = host_ens;
		_dens = device_ens;
	}

	//! Synchronize device ensemble with host ensemble
	void upload_ensemble() {
		_hens.copyTo( _dens );
	}
		
	//! Synchronize host ensemble with device ensemble
	void download_ensemble() {
		_dens.copyTo( _hens );
	}

	//! Grid dimentions for CUDA kernel launch
	virtual dim3 gridDim() = 0;
	//! Block dimentions for CUDA kernel launch
	virtual dim3 threadDim() = 0;
	//! Amount of shared memory required for CUDA kernel launch
	virtual int  shmemSize() = 0;
};
typedef shared_ptr<integrator> Pintegrator;


}

}
