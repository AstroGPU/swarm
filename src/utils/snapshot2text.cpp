#include "swarm/snapshot.hpp"
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
void parse_cmd(int argc, char* argv[], string& ifn, string& ofn){
	namespace po = boost::program_options;
	po::positional_options_description pos;
	po::options_description desc(std::string("Usage: ") + argv[0] + " \nOptions");

	desc.add_options()
		("input,i", po::value<std::string>(), "Input file")
		("output,o", po::value<std::string>(), "Output file")
		("help,h", "produce help message")
		;

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
	po::notify(vm);

	//// Respond to switches 
	//
	if ((vm.count("help")>0) || (vm.count("input") == 0) || (vm.count("output") == 0)) { std::cout << desc << "\n"; exit(1); }

	ifn = vm["input"].as<string>();
	ofn = vm["output"].as<string>();

}

int main(int argc,char * argv[]){
	string ifn, ofn;
	parse_cmd(argc,argv,ifn,ofn);
	swarm::defaultEnsemble ens = swarm::snapshot::load(ifn);
	swarm::snapshot::save_text(ens,ofn);
}
