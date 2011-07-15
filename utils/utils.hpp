/*************
 *  Author : Saleh Dindar
 *
 *
 */
#include "swarm.h"
#include "ensemble.hpp"

#define $__(x,line) (std::cerr << __FUNCTION__ << ":" << line << ": " << #x <<  " = " << (x) << std::endl)
#define DEBUG_OUTPUT(level,message) ( (DEBUG_LEVEL >= level) ? (std::cerr << __FUNCTION__ << ":" << __LINE__ << ": " << (message) << std::endl) : std::cerr )


void generate_ensemble(swarm::config& cfg, swarm::cpu_ensemble& ens)  ;
bool validate_configuration(swarm::config& cfg);
double find_max_energy_conservation_error(swarm::ensemble& ens, swarm::ensemble& reference_ensemble ) ;
swarm::hostEnsemble generate_ensemble(swarm::config& cfg)  ;
void outputConfigSummary(std::ostream& o,swarm::config& cfg);
