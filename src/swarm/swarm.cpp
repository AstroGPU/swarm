#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include "swarm.h"
#include "query.hpp"
#include "snapshot.hpp"
#include "stopwatch.h"
#include "gpu/device_settings.hpp"

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
	double destination_time = cfg.optional("destination_time", begin_time + 10 * M_PI );
	double interval = cfg.optional("interval", (destination_time-begin_time) ) ; 
	double logarithmic = cfg.optional("logarithmic", 0 ) ; 
	double allowed_deltaE =  cfg.optional("allowed_deltaE", 0.01 );

	if(destination_time < begin_time ) ERROR("Destination time should be larger than begin time");
	if(interval < 0) ERROR("Interval cannot be negative");
	if(interval < 0.001) ERROR("Interval too small");
	if(logarithmic != 0 && logarithmic <= 1) ERROR("logarithm base should be greater than 1");


	std::cout << "Time, Energy Conservation Factor (delta E/E), Active Systems" << std::endl;
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
                ensemble::range_t deltaE_range = energy_conservation_error_range(ens, initial_ens );

		int active_systems = ens.nsys() - number_of_disabled_systems( ens );
		std::cout << effective_time << ", " << deltaE_range.max << ", " << deltaE_range.median << ", " << active_systems << std::endl;
		
		if(deltaE_range.median > allowed_deltaE){
			INFO_OUTPUT(0, "At least half of systems are too unstable to conserve energy. dE/E: " << deltaE_range.median << std::endl );
			exit(1);
		}

	}

}

bool read_input_file(defaultEnsemble& ens, config& cfg){

	if(cfg.valid("input") ) {
		INFO_OUTPUT(1, "Loading from binary file " << cfg["input"]);
		ens = swarm::snapshot::load(cfg["input"]);	
		INFO_OUTPUT(1,", time = " << ens.time_ranges() << endl);
		return true;

	}else if(cfg.valid("text_input")) {
		INFO_OUTPUT(1, "Loading from text file " << cfg["text_input"]);
		ens = swarm::snapshot::load_text(cfg["text_input"]);	
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

	}else if(cfg.valid("text_output")) {
		INFO_OUTPUT(1, "Loading from text file " << cfg["text_output"]);
		ens = swarm::snapshot::load_text(cfg["text_output"]);	
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
		srand(time(NULL));
		initial_ens = generate_ensemble(cfg);
	}
}

bool save_ensemble(){ 
	// Save the ensemble
	if(cfg.valid("output")) {
		INFO_OUTPUT(1, "Saving as binary  " << cfg["output"]);
		INFO_OUTPUT(1, ", time = " << current_ens.time_ranges() << endl);
		swarm::snapshot::save(current_ens,cfg["output"]);	
		return true;

	}else if(cfg.valid("text_output")) {
		INFO_OUTPUT(1, "Saving as text  " << cfg["text_output"]);
		INFO_OUTPUT(1, ", time = " << current_ens.time_ranges() << endl);
		swarm::snapshot::save_text(current_ens,cfg["text_output"]);	
		return true;

	}else
		return false;
}

void init_cuda(){
	// Initialize Swarm
	swarm::init(cfg);
	//	print_device_information();
}

void prepare_integrator () {
	// Initialize Integrator
	DEBUG_OUTPUT(2, "Initializing integrator" );
	double begin_time = initial_ens.time_ranges().average;
	double destination_time = cfg.optional("destination_time", begin_time + 10 * M_PI );
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

		INFO_OUTPUT(1,"\tPosition difference: " << pos_diff  << endl
			 << "\tVelocity difference: " << vel_diff  << endl
			 << "\tTime     difference: " << time_diff << endl );

		if( !comparison || pos_diff > pos_threshold || vel_diff > vel_threshold || time_diff > time_threshold ){
			INFO_OUTPUT(0, "Test failed" << endl);
			exit(1);
		}else {
			INFO_OUTPUT(0, "Test success" << endl);
		}


	}else{
		ERROR("You should provide a test output file");
	}

	INFO_OUTPUT( 1, "Integration time: " << integration_time << " ms " << std::endl);
}

void benchmark_item(const string& param, const string& value) {
	//outputConfigSummary(std::cout,cfg);
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

/**
 * Data structure hold the values for a range
 * of parameters
 */
struct parameter_range {
	string parameter;
	vector<string> values;
	bool interpolate;
};

/// Because the signature of the function changed after boost version 1.41
/// we have to put this guard here
#if BOOST_VERSION  < 104200
	void raise_validation_error(const string& s){
		throw po::validation_error(s);
	}
#else
	void raise_validation_error(const string& s){
		throw po::validation_error(po::validation_error::invalid_option_value, s);
	}
#endif

/// Function called by Boost to parse the range parameter
void validate(boost::any& v, const std::vector<std::string>& values
	, parameter_range* , int){

	// Make sure no previous assignment to 'v' was made.
	po::validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	const std::string& s = po::validators::get_single_string(values);

	parameter_range p;
	static boost::regex assign("^(\\w+)=(.+)$");
	static boost::regex range("^([0-9\\.Ee]+)\\.\\.([0-9\\.Ee]+)$");
	static boost::regex rangeinc("^([0-9\\.Ee]+)\\.\\.([0-9\\.Ee]+)\\.\\.([0-9\\.Ee]+)$");
	static boost::regex list("^([a-zA-Z0-9\\.]+)(?:,([a-zA-Z0-9\\.]+))+$");
	static boost::regex item("([a-zA-Z0-9\\.]+)");
	boost::smatch match, items;
	if (boost::regex_match(s, match, assign)) {
		p.parameter = match[1];
		string value_range = match[2];
		boost::smatch range_match, list_match;
		if (boost::regex_match( value_range , range_match, rangeinc )) {
			p.values.push_back(range_match[1]);
			p.values.push_back(range_match[3]);
			p.values.push_back(range_match[2]);
			p.interpolate = true;
		}else if (boost::regex_match( value_range , range_match, range )) {
			p.values.push_back(range_match[1]);
			p.values.push_back(range_match[2]);
			p.interpolate = true;
		}else if (boost::regex_match( value_range , list_match, list )) {
			boost::sregex_iterator items( value_range.begin(), value_range.end() , item ), end;
			for(;items != end; items++)
				p.values.push_back((*items)[0]);
			p.interpolate = false;
		}else {
			raise_validation_error("Range is invalid");
		}
	}else
		raise_validation_error("Wrong parameter-value pair ");
	v = boost::any(p);
}

void benchmark(const parameter_range& pr){ 


	// Go through the loop
	std::cout << "Parameter, Value, Time, Energy Conservation Factor (delta E/E)"
				 ", Position Difference, Velocity Difference, Time difference "
			     ", Integration (ms), Integrator initialize (ms) \n";

	string param = pr.parameter;
	vector<string> values = pr.values;

	if(!(param == "input" || param == "nsys" || param == "nbod" || cfg.count("reinitialize")))
		load_generate_ensemble();
	outputConfigSummary(std::cout,cfg);
	if(verify_mode)
		reference_integration();

	if(pr.interpolate)  {
		// Numeric range
		double from = atof(values[0].c_str()), to = atof(values[1].c_str());
		double increment = (values.size() == 3) ? atof(values[2].c_str()) : 1;

		for(double i = from; i <= to ; i+= increment ){
			std::ostringstream stream;
			stream <<  i;
			string value = stream.str();

			cfg[param] = value;

			DEBUG_OUTPUT(1, "=========================================");
			benchmark_item(param,value);
		}
	}else {
		// Iterate over values
		for(vector<string>::const_iterator i = values.begin();  i != values.end(); i++){
			string value = *i;

			if(param == "config")
				cfg = config::load(value,base_cfg);
			else
				cfg[param] = *i;

			DEBUG_OUTPUT(1, "=========================================");
			benchmark_item(param,value);
		}
	}
}

void print_version(){
	bool amd64 = sizeof(void*) == 8;
	cout << "Swarm is running as " << (amd64 ? "64bit" : "32bit" ) << endl;
	exit(0);
}

void parse_commandline_and_config(int argc, char* argv[]){

	po::positional_options_description pos;
	pos.add("command", 1);
	pos.add("param_value_pair", -1);

	po::options_description desc("Usage:\n \tswarm [options] COMMAND [PARAMETER VALUE VALUE ...]\n\n"
			"Possible commands are: \n"
			"\tintegrate :  Integrate a [loaded|generated] ensemble\n"
			"\tbenchmark :  Compare outputs for different methods of integrations\n"
			"\tverify    :  Verify an integrator against a reference integrator\n"
			"\tquery     :  Query data from a log file\n"
			"\ttest      :  Test a configuration against input/output files\n"
			"\ttest-cpu  :  Test a configuration against input/output files without initializing GPU\n"
			"\tgenerate  :  Generate a new ensemble and save it to output file\n"
			"\tconvert   :  Read input file and write it to output file (converts to/from text)\n"
			"\nOptions"
		);


	po::options_description integrate("Integation Options");
	integrate.add_options()
		("destination_time,d", po::value<std::string>() , "Destination time to achieve")
		("logarithmic,l", po::value<std::string>() , "Produce times in logarithmic scale" )
		("interval,n", po::value<std::string>() , "Energy test intervals")
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("text_input,I", po::value<std::string>(), "Text Input file")
		("text_output,O", po::value<std::string>(), "Text Output file");

	po::options_description benchmark("Benchmark Options");
	benchmark.add_options()
		("range", po::value<parameter_range>(), "Parameter value range param=v1,v2,v3");

	po::options_description query("Query Options");
	query.add_options()
		("time,t", po::value<time_range_t>(), "range of times to query")
		("system,s", po::value<sys_range_t>(), "range of systems to query")
		("body,b", po::value<sys_range_t>(), "range of bodies to query")
		("keplerian,k", "output in Keplerian coordinates")
		("astrocentric", "output coordinates in astrocentric frame")
		("barycentric", "output coordinates in barycentric frame")
		("origin", "output coordinates in origin frame [default w/ Cartesian]")
		("jacobi", "output coordinates in Jacobi frame [default w/ Keplerian]")
		("logfile,f", po::value<std::string>(), "the log file to query");

	po::options_description positional("Positional Options");
	positional.add_options()
		("command" , po::value<string>() , "Swarm command: integrate, benchmark, verify, query ")
		("param_value_pair" , po::value< vector<string> >() , "Parameter value pairs: param=value")
		;

	po::options_description general("General Options");
	general.add_options()
		("cfg,c", po::value<std::string>(), "Integrator configuration file")
		("defaults", "Use default values for required configuration options")
		("help,h", "produce help message")
		("plugins,p", "list all of the plugins")
		("verbose,v", po::value<int>(), "Verbosity level (debug output) ")
		("quiet,q",  "Suppress all the notifications (for CSV generation)")
		("version,V",  "Print version message");
	
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
	if (vm.count("version")) print_version();
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
	const char* cmd_to_config[cmd_to_config_len] = { "input", "output", "text_input" , "text_output", "destination_time", "logarithmic", "interval" };
	for(int i = 0; i < cmd_to_config_len; i++)
		if(vm.count(cmd_to_config[i]))
			cfg[cmd_to_config[i]] = vm[cmd_to_config[i]].as<std::string>();

	if(vm.count("param_value_pair") > 0) {
		vector<string> pairs = vm["param_value_pair"].as< vector<string> >();

		for(vector<string>::const_iterator i = pairs.begin();  i != pairs.end(); i++){
			string::size_type pos     = i->find_first_of( "=" );
			if( string::npos != pos ){
				string parameter = i->substr( 0, pos );
				string value = i->substr( pos + 1 );
				cfg[parameter] = value;
				INFO_OUTPUT(3, parameter << " set to " << value << endl );
			}else{
				ERROR( ("Wrong parameter value pair : " + *i ).c_str() );
			}
		}
	}


}


int main(int argc, char* argv[]){

	// Command line and Config file
	parse_commandline_and_config(argc,argv);
//	outputConfigSummary(std::cout,cfg);
	base_cfg = cfg;
	command = argvars_map["command"].as< string >();


	// Set some variables
	pos_threshold = cfg.optional("pos_threshold", 1e-10);
	vel_threshold = cfg.optional("vel_threshold", 1e-10);
	time_threshold = cfg.optional("time_threshold", 1e-4);

	// Branch based on COMMAND
	if(command == "integrate"){
		init_cuda();
		run_integration();

	}else if(command == "test-cpu"){
//		init_cuda();  // removed so CPU tests would work, but then broke GPU tests when needed to set cuda device
		output_test();
	}else if(command == "test"){
		init_cuda();  
		output_test();
	}else if(command == "benchmark" || command == "verify") {
		init_cuda();
		if(command == "verify") 
			verify_mode = true;

		if(argvars_map.count("range") > 0) {
			parameter_range pr = argvars_map["range"].as<parameter_range>();
			benchmark(pr);
		}else
			cerr << "A parameter-value range is missing" << endl;

		if(command == "verify")  {
			cout << (verification_results ? "Verify success" : "Verify failed") << endl;
			return verification_results ? 0 : 1;
		}
	}


	else if(command == "query" ) {
		if (!argvars_map.count("logfile")) { cerr << "Name of input log file is missing \n"; return 1; }

		time_range_t T;
		sys_range_t sys, body_range;
		if (argvars_map.count("time")) { T = argvars_map["time"].as<time_range_t>(); }
		if (argvars_map.count("system")) { sys = argvars_map["system"].as<sys_range_t>(); }
		if (argvars_map.count("body")) { body_range = argvars_map["body"].as<sys_range_t>(); }
		if (argvars_map.count("keplerian")) { query::set_keplerian_output(); }
		if (argvars_map.count("origin")) { query::set_coordinate_system(query::origin); }
		if (argvars_map.count("astrocentric")) { query::set_coordinate_system(query::astrocentric); }
		if (argvars_map.count("barycentric")) { query::set_coordinate_system(query::barycentric); }
		if (argvars_map.count("jacobi")) { query::set_coordinate_system(query::jacobi); }


        std::cout.flush(); 
	cerr << "# Systems: " << sys << " Bodies: " << body_range << " Times: " << T << endl;
	if (argvars_map.count("keplerian")) 
	  { std::cerr << "# Output in Keplerian coordinates  "; }
	else { std::cerr << "# Output in Cartesian coordinates  "; }
#if 0
	switch(query::planets_coordinate_system)
	  {
		case astrocentric:
		  std::cerr << "(astrocentric)";
		  break;
		case barycentric:
		  std::cerr << "(barycentric)";
		  break;
		case jacobi:
		  std::cerr << "(jacobi)";
		  break;
		case origin:
		  std::cerr << "(origin)";
		  break;
	  }
#endif
	std::cerr << "\n";


		std::string datafile(argvars_map["logfile"].as<std::string>());

		query::execute(datafile, T, sys, body_range);
	}

	else if(command == "generate" ) {
		srand(time(NULL));
		INFO_OUTPUT(1, "Generating new ensemble:  " << cfg["nsys"] << ", " << cfg["nbod"] << endl);
		current_ens = generate_ensemble(cfg);
		save_ensemble();
	} 
	
	else if(command == "convert" ) {
		if( read_input_file(current_ens,cfg) && save_ensemble() ) {
			INFO_OUTPUT(1,"Converted!");
		}else{
			INFO_OUTPUT(1,"Convertion failed");
			exit(1);
		}
	} 
	
	else
		std::cerr << "Valid commands are: integrate, benchmark, verify, test, test-cpu, query, generate " << std::endl;

	return 0;
}
