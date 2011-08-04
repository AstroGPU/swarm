#include "plugin_manager.hpp"
#include <map>
#include <string>

#include <iostream>
using namespace std;

namespace swarm {

typedef map<string,plugin*> plugins_map_t;


plugins_map_t& plugins_map() {
	static plugins_map_t pmap;
	return pmap;
}

void add_plugin(plugin* p){
	plugins_map()[p->id()] = p;
}

plugin* get_plugin(const std::string& name){
	plugin* p = plugins_map()[name];
	if(p==0)
		throw plugin_not_found(name);
	return p;
}
void* instance_plugin(const std::string& name,const config& cfg){
	return get_plugin(name)->create(cfg);
}

vector<plugin*> get_all_plugins() {
	vector<plugin*> v;
	for(plugins_map_t::iterator i = plugins_map().begin(); i != plugins_map().end(); i++) {
		v.push_back(i->second);
	}
	return v;
}


plugin_help_message_t plugin_help_message;

ostream& operator << (ostream& out, const plugin_help_message_t&){
	vector<plugin*> ps = get_all_plugins();
	out << "Found " << ps.size() << " plugins " << endl;
	for(int i = 0; i < ps.size(); i++){
		out << ps[i]->id() << "\t" << ps[i]->description() << endl;
	}
}

}
