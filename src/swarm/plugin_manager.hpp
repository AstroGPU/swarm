#include "plugin.hpp"

#include <vector>
#include <ostream>
namespace swarm {

void add_plugin(plugin* p);
plugin* get_plugin(const std::string& name);
void* instance_plugin(const std::string& name, const config& cfg);
std::vector<plugin*> get_all_plugins();

struct plugin_help_message_t {};
extern plugin_help_message_t plugin_help_message;
std::ostream& operator << (std::ostream&, const plugin_help_message_t&);


struct plugin_not_found : std::exception {
	std::string _name;
	plugin_not_found(std::string name) : _name(name) {}
	virtual ~plugin_not_found() throw(){}

	virtual const char * what() const throw() { 
		return ("Plugin " + _name + " was not found ").c_str(); 
	}
};

}
