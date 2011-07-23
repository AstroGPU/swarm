/*! \file user.h
 *  \brief #defines that advanced users might want to change
 *
*/ 

// currently setup for swarm_adap.cpp
#define OUTPUT_DIRECTORY "adap_output"
#define LARGE_NUMBER 1e32
#define SMALL_NUMBER 1e-32

//uncomment for N bodies > 10
//#define __LARGE_N_SYSTEM__

#ifdef __LARGE_N_SYSTEM__
// N_BODY should be larger than 10. 12 is just an example 
#define N_BODY 12
#endif
