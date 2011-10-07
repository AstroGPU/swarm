#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "swarm.h"
#include "query.hpp"
#include "snapshot.hpp"
#include "stopwatch.h"

int DEBUG_LEVEL  = 0;

using namespace swarm;
using namespace std;
namespace po = boost::program_options;
using boost::bind;
using boost::shared_ptr;



// Runtime variables
string command;
config cfg;
config base_cfg;
defaultEnsemble initial_ens;
defaultEnsemble current_ens;
defaultEnsemble reference_ens;
po::variables_map argvars_map;
Pintegrator integ;

void stability_test() {

	defaultEnsemble &ens = current_ens;

	double begin_time = ens.time_ranges().median;
	double destination_time = cfg.optional("destination time", begin_time + 10 * M_PI );
	double interval = cfg.optional("interval", (destination_time-begin_time) ) ; 
	double logarithmic = cfg.optional("logarithmic", 0 ) ; 

	if(destination_time < begin_time ) ERROR("Destination time should be larger than begin time");
	if(interval < 0) ERROR("Interval cannot be negative");
	if(interval < 0.001) ERROR("Interval too small");
	if(logarithmic != 0 && logarithmic <= 1) ERROR("logarithm base should be greater than 1");


	std::cout << "Time, Energy Conservation Factor (delta E/E)" << std::endl;
	for(double time = begin_time; time < destination_time; ) {

		if(logarithmic > 1)
			time = (time > 0) ? time * logarithmic : interval;
		else
			time += interval;

		double effective_time = min(time,destination_time);
		integ->set_destination_time ( effective_time );

		DEBUG_OUTPUT(2, "Integrator ensemble" );
		integ->integrate();

		SYNC;
		DEBUG_OUTPUT(2, "Check energy conservation" );
		double max_deltaE = find_max_energy_conservation_error(ens, initial_ens );
		std::cout << effective_time << ", " << max_deltaE << std::endl;

	}

}

bool read_input_file(defaultEnsemble& ens, config& cfg){

	if(cfg.valid("input") ) {
		INFO_OUTPUT(1, "Loading from binary file " << cfg["input"]);
		ens = swarm::snapshot::load(cfg["input"]);	
		INFO_OUTPUT(1,", time = " << ens.time_ranges() << endl);
		return true;

	}else if(cfg.valid("text input")) {
		INFO_OUTPUT(1, "Loading from text file " << cfg["text input"]);
		ens = swarm::snapshot::load_text(cfg["text input"]);	
		INFO_OUTPUT(1,", time = " << ens.time_ranges() << endl);
		return true;

	}else
		return false;
		
}

bool read_output_file(defaultEnsemble& ens, config& cfg){

	if(cfg.valid("output") ) {
		INFO_OUTPUT(1, "Loading from binary file " << cfg["output"]);
		ens = swarm::snapshot::load(cfg["output"]);	
		INFO_OUTPUT(1,", time = " << ens.time_ranges() << endl);
		return true;

	}else if(cfg.valid("text output")) {
		INFO_OUTPUT(1, "Loading from text file " << cfg["text output"]);
		ens = swarm::snapshot::load_text(cfg["text output"]);	
		INFO_OUTPUT(1,", time = " << ens.time_ranges() << endl);
		return true;

	}else
		return false;
		
}

void load_generate_ensemble(){
	// Load/Generate the ensemble
	if(read_input_file(initial_ens,cfg)) {
		
	}else{
		INFO_OUTPUT(1, "Generating new ensemble:  " << cfg["nsys"] << ", " << cfg["nbod"] << endl);
		initial_ens = generate_ensemble(cfg);
	}
}

void save_ensemble(){ 
	// Save the ensemble
	if(cfg.valid("output")) {
		INFO_OUTPUT(1, "Saving as binary  " << cfg["output"]);
		INFO_OUTPUT(1, ", time = " << current_ens.time_ranges() << endl);
		swarm::snapshot::save(current_ens,cfg["output"]);	

	}else if(cfg.valid("text output")) {
		INFO_OUTPUT(1, "Saving as text  " << cfg["text output"]);
		INFO_OUTPUT(1, ", time = " << current_ens.time_ranges() << endl);
		swarm::snapshot::save_text(current_ens,cfg["text output"]);	
	}
}

void prepare_integrator () {
	// Initialize Integrator
	DEBUG_OUTPUT(2, "Initializing integrator" );
	double begin_time = initial_ens.time_ranges().average;
	double destination_time = cfg.optional("destination time", begin_time + 10 * M_PI );
	integ = integrator::create(cfg);
	integ->set_ensemble(current_ens);
	integ->set_destination_time ( destination_time );
	SYNC;
}

void generic_integrate () {
	DEBUG_OUTPUT(2, "Integrator ensemble" );
	integ->integrate();
}

void run_integration(){
	if(!validate_configuration(cfg) ) ERROR( "Invalid configuration" );

	load_generate_ensemble();

	DEBUG_OUTPUT(2, "Make a copy of ensemble for energy conservation test" );
	current_ens = initial_ens.clone();

	prepare_integrator();

	double integration_time = watch_time ( cfg.valid("interval") ? stability_test : generic_integrate );

	save_ensemble();

	INFO_OUTPUT( 1, "Integration time: " << integration_time << " ms " << std::endl);
}


void reference_integration() {
	DEBUG_OUTPUT(2, "Make a copy of ensemble for reference ensemble" );
	reference_ens = initial_ens.clone();
	DEBUG_OUTPUT(1, "Reference integration" );
	prepare_integrator();
	integ->set_ensemble( reference_ens );
	integ->integrate();
}



bool verification_results = true, verify_mode = false;
double pos_threshold = 0, vel_threshold = 0, time_threshold = 0;

void fail_verify() {
	if(verify_mode){
		INFO_OUTPUT(0, "Verify failed" << endl);
		exit(1);
	}else{
		verification_results = false;
	}
}

void output_test() {
	if(!validate_configuration(cfg) ) ERROR( "Invalid configuration" );

	if(!read_input_file(initial_ens, cfg) ) {
		ERROR("you should have a tested input file");
	}

	DEBUG_OUTPUT(2, "Make a copy of ensemble for energy conservation test" );
	current_ens = initial_ens.clone();

	prepare_integrator();

	double integration_time = watch_time ( cfg.valid("interval") ? stability_test : generic_integrate );

	if(read_output_file(reference_ens,cfg)){

		// Compare with reneference ensemble for integrator verification
		double pos_diff = 0, vel_diff = 0, time_diff = 0;
		bool comparison =  compare_ensembles( current_ens, reference_ens , pos_diff, vel_diff, time_diff );

		if( !comparison || pos_diff > pos_threshold || vel_diff > vel_threshold || time_diff > time_threshold ){
			INFO_OUTPUT(0, "Test failed" << endl);
		}else {
			INFO_OUTPUT(0, "Test success" << endl);
		}

		INFO_OUTPUT(1,"\tPosition difference: " << pos_diff  << endl
			 << "\tVelocity difference: " << vel_diff  << endl
			 << "\tTime     difference: " << time_diff << endl );

	}else{
		ERROR("You should provide a test output file");
	}

	INFO_OUTPUT( 1, "Integration time: " << integration_time << " ms " << std::endl);
}

void benchmark_item(const string& param, const string& value) {
	if(!validate_configuration(cfg) ) ERROR( "Invalid configuration" );

	if(param == "input" || param == "nsys" || param == "nbod" || cfg.count("reinitialize"))
		load_generate_ensemble();

	DEBUG_OUTPUT(2, "Make a copy of ensemble for energy conservation test" );
	current_ens = initial_ens.clone();

	double init_time = watch_time( prepare_integrator );

	double integration_time = watch_time( generic_integrate );

	double max_deltaE = find_max_energy_conservation_error(current_ens, initial_ens );

	// Compare with reneference ensemble for integrator verification
	double pos_diff = 0, vel_diff = 0, time_diff = 0;
	bool comparison =  compare_ensembles( current_ens, reference_ens , pos_diff, vel_diff, time_diff );

	/// CSV output for use in spreadsheet software 
	std::cout << param << ", "
	          << value << ",  "   
			  << current_ens.time_ranges() << ",  "
	          << max_deltaE << ",    " 
			  << pos_diff << ",    "
			  << vel_diff << ",   "
			  << time_diff << ",   "
	          << integration_time << ",    "
	          << init_time 
	          << std::endl;

	if( !comparison || pos_diff > pos_threshold || vel_diff > vel_threshold || time_diff > time_threshold ){
		fail_verify();
	}

}

void benchmark(){ 
	// Reference settings
	if(!validate_configuration(cfg) ) ERROR( "Invalid configuration" );
	load_generate_ensemble();
	reference_integration();


	// Go through the loop
	std::cout << "Parameter, Value, Time, Energy Conservation Factor (delta E/E)"
				 ", Position Difference, Velocity Difference, Time difference "
			     ", Integration (ms), Integrator initialize (ms) \n";

	po::variables_map &vm = argvars_map;
	if((vm.count("parameter") > 0) && (vm.count("value") > 0)) {
		
		// Iterate over values
		string param = vm["parameter"].as< vector<string> >().front();
		vector<string> values = vm["value"].as< vector<string> >();

		for(vector<string>::const_iterator i = values.begin();  i != values.end(); i++){
			string value = *i;

			if(param == "config")
				cfg = config::load(value,base_cfg);
			else
				cfg[param] = *i;

			DEBUG_OUTPUT(1, "=========================================");
			benchmark_item(param,value);
		}

	} else if ((vm.count("parameter") > 0) &&(vm.count("from") > 0) && (vm.count("to") > 0)){

		// Numeric range
		int from = vm["from"].as<int>(), to = vm["to"].as<int>();
		int increment = vm.count("inc") > 0 ? vm["inc"].as<int>() : 1;
		string param = vm["parameter"].as< vector<string> >().front();

		for(int i = from; i <= to ; i+= increment ){
			std::ostringstream stream;
			stream <<  i;
			string value = stream.str();

			cfg[param] = value;

			DEBUG_OUTPUT(1, "=========================================");
			benchmark_item(param,value);
		}

	} else {
		std::cerr << "Invalid combination of parameters and values " << std::endl;
	}
}

void parse_commandline_and_config(int argc, char* argv[]){

	po::positional_options_description pos;
	pos.add("command", 1);
	pos.add("parameter", 1);
	pos.add("value", -1);

	po::options_description desc("Usage:\n \tswarm [options] COMMAND [PARAMETER VALUE VALUE ...]\n\n"
			"Possible commands are: \n"
			"\tintegrate :  Integrate a [loaded|generated] ensemble\n"
			"\tbenchmark :  Compare outputs for different methods of integrations\n"
			"\tverify    :  Verify an integrator against a reference integrator\n"
			"\tquery     :  Query data from a log file\n"
			"\ttest      :  Test a configuration against input/output files\n"
			"\tgenerate  :  Generate a new ensemble and save it to output file\n"
			"\tconvert   :  Read input file and write it to output file (converts to/from text)\n"
			"\nOptions"
		);


	po::options_description integrate("Integation Options");
	integrate.add_options()
		("destination time,d", po::value<std::string>() , "Destination time to achieve")
		("logarithmic,l", po::value<std::string>() , "Produce times in logarithmic scale" )
		("interval,n", po::value<std::string>() , "Energy test intervals")
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("text input,I", po::value<std::string>(), "Text Input file")
		("text output,O", po::value<std::string>(), "Text Output file");

	po::options_description benchmark("Benchmark Options");
	benchmark.add_options()
		("from", po::value<int>() , "from integer value")
		("to", po::value<int>() , "to integer value")
		("inc", po::value<int>() , "increment integer value");

	po::options_description query("Query Options");
	query.add_options()
		("time,t", po::value<time_range_t>(), "range of times to query")
		("system,s", po::value<sys_range_t>(), "range of systems to query")
		("logfile,f", po::value<std::string>(), "the log file to query");

	po::options_description positional("Positional Options");
	positional.add_options()
		("command" , po::value<string>() , "Swarm command: integrate, benchmark, verify, query ")
		("parameter" , po::value<vector<string> >() , "Parameteres to benchmark: config, integrator, nsys, nbod, blocksize, ... ")
		("value" , po::value<vector<string> >() , "Values to iterate over ");

	po::options_description general("General Options");
	general.add_options()
		("cfg,c", po::value<std::string>(), "Integrator configuration file")
		("defaults", "Use default values for required configuration options")
		("help,h", "produce help message")
		("plugins,p", "list all of the plugins")
		("verbose,v", po::value<int>(), "Verbosity level (debug output) ")
		("quiet,q",  "Suppress all the notifications (for CSV generation)");
	
	desc.add(general).add(integrate).add(benchmark).add(query);

	po::options_description all;
	all.add(desc).add(positional);
	


	po::variables_map &vm = argvars_map;
	po::store(po::command_line_parser(argc, argv).
			options(all).positional(pos).run(), vm);
	po::notify(vm);

	//// Respond to switches 
	//
	if (vm.count("help")) { std::cout << desc << "\n"; exit(1); }
	if (vm.count("verbose") ) DEBUG_LEVEL = vm["verbose"].as<int>();
	if (vm.count("quiet") ) DEBUG_LEVEL = -1;

	if (vm.count("plugins")) {
		cout << swarm::plugin::help;
		exit(0);
	}

	if(vm.count("defaults")){
		cfg = default_config();
	}
	if(vm.count("cfg")){
		std::string icfgfn =  vm["cfg"].as<std::string>();
		cfg = config::load(icfgfn, cfg);
	}

	if(vm.count("command") == 0){
		std::cout << desc << "\n"; exit(1); 
	}

	const int cmd_to_config_len = 7;
	const char* cmd_to_config[cmd_to_config_len] = { "input", "output", "text input" , "text output", "destination time", "logarithmic", "interval" };
	for(int i = 0; i < cmd_to_config_len; i++)
		if(vm.count(cmd_to_config[i]))
			cfg[cmd_to_config[i]] = vm[cmd_to_config[i]].as<std::string>();


}


int main(int argc, char* argv[]){

	// Command line and Config file
	parse_commandline_and_config(argc,argv);
//	outputConfigSummary(std::cout,cfg);
	base_cfg = cfg;
	command = argvars_map["command"].as< string >();

	// Initialize Swarm
	swarm::init(cfg);
	srand(time(NULL));

	// Set some variables
	pos_threshold = cfg.optional("pos threshold", 1e-10);
	vel_threshold = cfg.optional("vel threshold", 1e-10);
	time_threshold = cfg.optional("time threshold", 1e-4);

	// Branch based on COMMAND
	if(command == "integrate")
		run_integration();

	else if(command == "test")
		output_test();

	else if(command == "benchmark") {
		benchmark();
	}

	else if(command == "verify" ) {
		verify_mode = true;
		benchmark();
		cout << (verification_results ? "Verify success" : "Verify failed") << endl;
		return verification_results ? 0 : 1;
	}

	else if(command == "query" ) {
		if (!argvars_map.count("logfile")) { cerr << "Name of input log file is missing \n"; return 1; }

		time_range_t T;
		sys_range_t sys;
		if (argvars_map.count("time")) { T = argvars_map["time"].as<time_range_t>(); }
		if (argvars_map.count("system")) { sys = argvars_map["system"].as<sys_range_t>(); }

        std::cout.flush(); 
		cerr << "Systems " << sys << " in time range " << T << endl;

		std::string datafile(argvars_map["logfile"].as<std::string>());

		query::execute(datafile, T, sys);
	}

	else if(command == "generate" ) {
		INFO_OUTPUT(1, "Generating new ensemble:  " << cfg["nsys"] << ", " << cfg["nbod"] << endl);
		current_ens = generate_ensemble(cfg);
		save_ensemble();
	} 
	
	else if(command == "convert" ) {
		if( read_input_file(current_ens,cfg) ) {
			save_ensemble();
			INFO_OUTPUT(1,"Converted!");
		}else{
			INFO_OUTPUT(1,"Convertion failed");
		}
	} 
	
	else
		std::cerr << "Valid commands are: integrate, benchmark, verify " << std::endl;

	return 0;
}
