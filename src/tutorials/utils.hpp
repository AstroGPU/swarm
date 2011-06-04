/*************
 *  Author : Saleh Dindar
 *
 *
 */
#include "swarm.h"
#include "hp/datatypes.hpp"

#define $_(x) (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": " << #x <<  " = " << (x) << std::endl)
//#define $$ (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": @ here " << std::endl)
#define $__(x,line) (std::cerr << __FUNCTION__ << ":" << line << ": " << #x <<  " = " << (x) << std::endl)
#define DEBUG_OUTPUT(level,message) ( (DEBUG_LEVEL >= level) ? (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": " << (message) << std::endl) : std::cerr )


void generate_ensemble(swarm::config& cfg, swarm::cpu_ensemble& ens)  ;
double find_max_energy_conservation_error(swarm::cpu_ensemble& ens, swarm::cpu_ensemble& reference_ensemble ) ;
bool validate_configuration(swarm::config& cfg);
double find_max_energy_conservation_error(swarm::hp::ensemble& ens, swarm::hp::ensemble& reference_ensemble ) ;
swarm::hp::ensemble generate_ensemble(swarm::config& cfg)  ;
