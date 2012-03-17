#include "swarm/swarm.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>


swarm::gpu::Pintegrator create_gpu_integrator(const swarm::config& cfg){
	return boost::dynamic_pointer_cast<swarm::gpu::integrator>(swarm::integrator::create(cfg));
}

#define PROPERTY_ACCESSOR(CLASS,PROPERTY,TYPE)       \
	void set_##PROPERTY( CLASS &r, const TYPE & t)   \
	{ r.PROPERTY() = t;  	}                        \
	TYPE get_##PROPERTY( CLASS &r)                   \
	{ return r.PROPERTY(); }

PROPERTY_ACCESSOR(swarm::ensemble::SystemRef, time, double )
PROPERTY_ACCESSOR(swarm::ensemble::SystemRef, id  , int    )

PROPERTY_ACCESSOR(swarm::ensemble::Body     , mass, double )
PROPERTY_ACCESSOR(swarm::ensemble::Body::Component, pos, double )
PROPERTY_ACCESSOR(swarm::ensemble::Body::Component, vel, double )

/* Swarm as a python module
 *
 * The reason for lib prefix is because CMake automatically
 * adds lib prefix to the name of the target
 *
 */
BOOST_PYTHON_MODULE(libswarmng_ext) {
	using namespace boost::python;
	using boost::noncopyable;
	using namespace swarm;

	def("init", swarm::init );
	def("generate_ensemble", generate_ensemble );
	def("sync", cudaThreadSynchronize );

	class_<config>("Config")
		.def( map_indexing_suite< config >() )
		;

	class_<ensemble::Sys::attributes_t, noncopyable >( "Sys.Attributes", no_init )
		.def("__getitem__", &ensemble::Sys::attributes_t::getitem )
		.def("__setitem__", &ensemble::Sys::attributes_t::setitem )
		;

	class_<ensemble::Body::attributes_t, noncopyable >( "Body.Attributes", no_init )
		.def("__getitem__", &ensemble::Body::attributes_t::getitem )
		.def("__setitem__", &ensemble::Body::attributes_t::setitem )
		;

	class_<ensemble::Body::Component >("Body.Components", no_init)
		.add_property("pos", &get_pos, &set_pos )
		.add_property("vel", &get_vel, &set_vel )
		;

	/* All the essential properties are covered
	 * However, we can make some convenience access
	 * methods like x,y,z,vy,vx,vz
	 * or get positions and velocities as a python list
	 * or set the body from a python list, etc.
	 *
	 *
	 */
	class_<ensemble::Body >("Body", no_init )
		.add_property("mass", &get_mass, &set_mass)
		.def("radius", &ensemble::Body::radius )
		.def("speed", &ensemble::Body::speed )
		.def("attributes", &ensemble::Body::attributes, return_value_policy<reference_existing_object>() )
		.def("__getitem__", &ensemble::Body::getitem, return_value_policy<reference_existing_object>() )
		//.add_property("pos", &get_pos_list, &set_pos_list ) Not implemented
		//.add_property("vel", &get_vel_list, &set_vel_list ) Not implemented
		;

	/*
	 *  All important methods are covered for this class
	 */
	class_<ensemble::SystemRef  >("SystemRef", no_init )
		.add_property("time", &get_time, &set_time)
		.add_property("id", &get_id, &set_id)
		.add_property("total_energy", &ensemble::SystemRef::total_energy)
		.def("is_active", &ensemble::SystemRef::is_active )
		.def("is_inactive", &ensemble::SystemRef::is_inactive )
		.def("is_enabled", &ensemble::SystemRef::is_enabled )
		.def("set_active", &ensemble::SystemRef::set_active )
		.def("set_inactive", &ensemble::SystemRef::set_inactive )
		.def("set_disabled", &ensemble::SystemRef::set_disabled )
		.def("__getitem__", &ensemble::SystemRef::getitem, return_value_policy<reference_existing_object>() )
		.def("attributes", &ensemble::SystemRef::attributes, return_value_policy<reference_existing_object>() )

		;

	class_<ensemble, noncopyable >("Ensemble", no_init )
		.def("__length__", &ensemble::nsys , return_value_policy<copy_const_reference>() )
		.add_property("nsys", make_function(&ensemble::nsys , return_value_policy<copy_const_reference>() ) )
		.add_property("nbod", make_function(&ensemble::nbod , return_value_policy<copy_const_reference>() ) )
		.def("__getitem__", &ensemble::getitem )
		;
	class_<defaultEnsemble, bases<ensemble> >("DefaultEnsemble")
		.def( "clone", &defaultEnsemble::clone )
		;
	

	class_<integrator, Pintegrator, noncopyable >("Integrator", no_init )
		.def("create",&integrator::create)
		.staticmethod("create")
		.def("integrate", &integrator::integrate )
		.add_property("ensemble", make_function(&integrator::get_ensemble, return_value_policy<reference_existing_object>() ), &integrator::set_ensemble)
		.add_property("destination_time", &integrator::get_destination_time, &integrator::set_destination_time)
		;

	void (gpu::integrator::*gpu_set_ensemble)(defaultEnsemble&) = &gpu::integrator::set_ensemble;

	class_<gpu::integrator, bases<integrator> , gpu::Pintegrator, noncopyable>("GpuIntegrator", no_init)
		.def("create", &create_gpu_integrator )
		.staticmethod("create")
		.def("integrate", &gpu::integrator::integrate )
		.def("core_integrate", &gpu::integrator::core_integrate )
		.def("download_ensemble", &gpu::integrator::download_ensemble )
		.def("upload_ensemble", &gpu::integrator::upload_ensemble )
		.add_property("ensemble", make_function(&integrator::get_ensemble, return_value_policy<reference_existing_object>() ), gpu_set_ensemble)
		;

	def("find_max_energy_conservation_error", find_max_energy_conservation_error );





}
