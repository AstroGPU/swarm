#include "swarm/swarm.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

const char*greet(){
	return "Hello World!";
}

std::string join_strings(std::vector<std::string> a){
	std::string ret = "";
	for(int i = 0; i < a.size(); i++){
		ret += a[i] + ", ";
	}
	return ret;
}

swarm::gpu::Pintegrator create_gpu_integrator(const swarm::config& cfg){
	return boost::dynamic_pointer_cast<swarm::gpu::integrator>(swarm::integrator::create(cfg));
}

BOOST_PYTHON_MODULE(libswarmng_ext) {
	using namespace boost::python;
	using namespace swarm;

	def("init", swarm::init );
	def("generate_ensemble", generate_ensemble );
	def("sync", cudaThreadSynchronize );

	class_<config>("Config")
		.def( map_indexing_suite< config >() )
		;

	class_<ensemble, boost::noncopyable >("Ensemble", no_init )
		;
	class_<defaultEnsemble, bases<ensemble> >("DefaultEnsemble")
		.def( "clone", &defaultEnsemble::clone )
		;
	

	class_<integrator, Pintegrator, boost::noncopyable >("Integrator", no_init )
		.def("create",&integrator::create)
		.staticmethod("create")
		.def("integrate", &integrator::integrate )
		.add_property("ensemble", make_function(&integrator::get_ensemble, return_value_policy<reference_existing_object>() ), &integrator::set_ensemble)
		.add_property("destination_time", &integrator::get_destination_time, &integrator::set_destination_time)
		;

	void (gpu::integrator::*gpu_set_ensemble)(defaultEnsemble&) = &gpu::integrator::set_ensemble;

	class_<gpu::integrator, bases<integrator> , gpu::Pintegrator, boost::noncopyable>("GpuIntegrator", no_init)
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
