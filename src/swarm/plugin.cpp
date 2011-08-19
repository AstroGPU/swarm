#include "common.hpp"
#include "plugin.hpp"

using namespace std;

namespace swarm {

typedef map<string,plugin*> plugins_map_t;


plugins_map_t& plugins_map() {
	static plugins_map_t pmap;
	return pmap;
}

void plugin::add(plugin* p){
	plugins_map()[p->id()] = p;
}

plugin* plugin::find(const std::string& name){
	plugin* p = plugins_map()[name];
	if(p==0)
		throw plugin_not_found(name);
	return p;
}
void* plugin::instance(const std::string& name,const config& cfg){
	return plugin::find(name)->create(cfg);
}

vector<plugin*> plugin::all() {
	vector<plugin*> v;
	for(plugins_map_t::iterator i = plugins_map().begin(); i != plugins_map().end(); i++) {
		v.push_back(i->second);
	}
	return v;
}


plugin::help_t plugin::help;

ostream& operator << (ostream& out, const plugin::help_t&){
	vector<plugin*> ps = plugin::all();
	out << "Found " << ps.size() << " plugins " << endl;
	for(int i = 0; i < ps.size(); i++){
		out << ps[i]->id() << "\t" << ps[i]->description() << endl;
	}
	return out;
}

}
