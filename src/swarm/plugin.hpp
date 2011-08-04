#include "config.hpp"

namespace swarm {

struct plugin {
	virtual void* create(const config& cfg) = 0;
	virtual std::string description() { return ""; }
	virtual std::string id() = 0;
};


template<class T> 
struct basic_plugin : public plugin {

	std::string _id, _description;
	basic_plugin(const std::string& i, const std::string& d = std::string() ) 
		: _id(i),_description(d) {}

	virtual void* create(const config& cfg) { return new T(cfg); }
	virtual std::string id() { return _id; }
	virtual std::string description() { return _description; }
};

void add_plugin(plugin*);

template<class T> 
struct plugin_initializer {
	T t;
	plugin_initializer() {
		add_plugin(&t);
	}
};

template<class T>
struct basic_plugin_initializer {
	basic_plugin<T> t;

	basic_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: t(id,description) {
		add_plugin(&t);
	}
};


template<class T>
struct integrator_plugin_initializer : public basic_plugin_initializer<T> {
	integrator_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: basic_plugin_initializer<T>("integrator_" + id,description) {}
};

template<class T>
struct writer_plugin_initializer : public basic_plugin_initializer<T> {
	writer_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: basic_plugin_initializer<T>("writer_" + id,description) {}
};

}
