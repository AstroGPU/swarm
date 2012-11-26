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

/*! \file plugin.hpp
 *   \brief Defines \ref swarm::plugin class and \ref swarm::basic_plugin class,
 *          implements the interface and utility functions for the plugin development 
 *          and plugin management system
 *
 *   Containts:
 *    - Interface classes for creating new plugins
 *    - Helper classes to add new plugins
 *    - Functions to instance plugins
 *
*/
#pragma once
#include "types/config.hpp"

namespace swarm {

/*! Abstract interface class for all plugins.
 *  
 *  Derived intstances should implement \ref create method
 *  for instancing the main class of the plugin package. \ref id
 *  is used to identify the package for the plugin look-up table.
 * 
 *  It is difficult to use the abstract class, one should use the helper
 * classes @ref integrator_plugin_initializer, @ref writer_plugin_initializer .
 */
struct plugin {
	//! Instance the plugin and return the main class of the package
	virtual void* create(const config& cfg) = 0;
	//! Human readable text to describe the purpose of the plugin
	virtual std::string description() { return ""; }
	//! Unique identifier for plugin, should use the correct prefix.
	//! For integrators use: "integrator_" prefix
	//! For writers use: "writer_" prefix
	virtual std::string id() = 0;

	//! Register a new plugin
	static void    add(plugin* p);
	//! Find a plugin by id()
	static plugin* find(const std::string& id);
	//! Find and instance the pluging with specified configuration. 
	static void*   instance(const std::string& name, const config& cfg);
	//! List of all registered plugins
	static std::vector<plugin*> all();

	struct help_t{};
	/*! Plugin the list of plugins in a nice format
	 *   To use:
	 *   \code
	 *     std::cout << plugin::help << std::endl;
	 * \endcode
	 */
	static help_t help;

};

/*! Template class for easy plugin development.
 *  
 *  class T should have a constructor that takes config as parameter.
 *
 *  Making plugins that make an instance of a specific class
 *  is a pretty common task. To define the plugin class for
 *  your plugin package:
 *   - Define one class that does all the work in your plugin. e.g. MySystem
 *   - Instantiate this template as follows:
 *  basic_plugin<MySystem>( "MySystem" , "This is my plugin that I use" )
 *  
 */
template<class T> 
struct basic_plugin : public plugin {

	std::string _id, _description;
	//! Specify id and description in the constructor
	basic_plugin(const std::string& i, const std::string& d = std::string() ) 
		: _id(i),_description(d) {}

	//! Instance the class T and return it
	virtual void* create(const config& cfg) { return new T(cfg); }
	//! id as specified in constructor
	virtual std::string id() { return _id; }
	//! description as specified in constructor
	virtual std::string description() { return _description; }
};

/*! Template class to add your plugin to the list of plugins.
 *
 * If your plugin class (derived from \ref plugin) is named MyPlugin,
 * to add it to swarm library you need to add this line to the global
 * scope of your source code:
 * \code
 *   plugin_initializer< MyPlugin > my_plugin;
 * \endcode
 *
 * For easy to use version look at \ref basic_plugin_initializer.
 */
template<class T> 
struct plugin_initializer {
	T t;
	plugin_initializer() {
		plugin::add(&t);
	}
};


/*! Template class for easy plugin development and management.
 *
 * Shorthand for plugin_initializer< basic_plugin < T > with some extra feautures.
 *
 * If you have created a plugin with main class as MyPlugin, and you want
 * to name the plugin as "integrator_MyPlugin", put this line in the global scope
 * of your source code:
 * \code
 *   basic_plugin_initializer< MyPlugin > my_plugin_initialize 
 *    ("integrator_MyPlugin", "This is my test plugin" );
 * \endcode
 */
template<class T>
struct basic_plugin_initializer {
	basic_plugin<T> t;

	//! Register the a basic_plugin<T> with specified id and description.
	basic_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: t(id,description) {
			plugin::add(&t);
	}
};


/*! Template to add a new integrator to swarm.
 *  
 *  Example: If you have developed a new integrator class "advanced_example",
 *  put this line in the global scope to add it to the swarmng library.
 *  The second string is the human readable briefing about your integrator
 *  \code
 *   integrator_plugin_initializer< advanced_example > advanced_example_initializer
 *   	("advanced_example" , "This is my advanced integrators that I implemented" );
 * \endcode
 */
template<class T>
struct integrator_plugin_initializer : public basic_plugin_initializer<T> {
	//! Register an integrator plugin with specified id and description. 
	//! Note: automatically adds the prefix "integrator_"
	integrator_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: basic_plugin_initializer<T>("integrator_" + id,description) {}
};

/*! Template to add a new writer plugin to swarm.
 *  
 *  This is rarely used. But if you have another output method
 *  for swarm log files. you can create a plugin and add it
 *  like this ( put in global scope ):
 *  \code
 *   writer_plugin_initializer< alternate_output > alternate_output_initializer
 *   	("alternate_output" , "This is my alternative output writer that I implemented" );
 * \endcode
 *  
 */
template<class T>
struct writer_plugin_initializer : public basic_plugin_initializer<T> {
	writer_plugin_initializer(const std::string& id
			, const std::string& description = std::string() ) 
		: basic_plugin_initializer<T>("writer_" + id,description) {}
};

/*! Thrown when a nonexisting plugin is requested.
 *  The name of the missing plugin is is included
 */
struct plugin_not_found : std::exception {
	//! name of the missing pluging
	std::string _name;
	plugin_not_found(std::string name) : _name(name) {}
	virtual ~plugin_not_found() throw(){}

	virtual const char * what() const throw() { 
		return ("Plugin \"" + _name + "\" was not found ").c_str(); 
	}
};

/*! Helper function to print out the plugin help message.
 *   Prints out a list of all registered plugins in swarmng library.
 *   To use:
 *   \code
 *     std::cout << plugin::help << std::endl;
 * \endcode
 *
 */
std::ostream& operator << (std::ostream&, const plugin::help_t&);

}
