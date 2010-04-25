README for Swarm-NG
===================
Swarm-NG Development Team
v0.1, April 2010: Initial release + 
Swarm-NG v0.1 (c) Swarm-NG Development Team 2010

Introduction
------------

.What is Swarm-NG?
Swarm-NG is a collection of code which to help scientists and engineers harness the power of GPUs.
In the early releases, Swarm-NG will focus on the integration of an ensemble of N-body systems evolving under Newtonian gravity.
There are existing libraries to perform the force calculation for large-N systems on GPUs.  
Swarm-NG does not replicate these, but rather focuses on integrating an ensemble of many systems where N is small. 
This is of partiular interest for astronomers who study the chaotic evolution of planetary systems.
In the long term, we hope Swarm-NG will allow for the efficient parallel integration of user-defined systems of ordinary differential equations.  

.What is included in {revdate} release of Swarm-NG {revnumber}?
This initial release contains GPU-based integrators for few-body systems, along with a draft of the public interface for setting initial conditions, calling the GPU-based integration algorithms and accessing the output data.  Small example programs are provided to demonstrate functionality and how to program with Swarm-NG.  

.Quick Start Tips
Read link:doc/build_system.html[doc/build_system.html] for how to build the library and example programs. + 
Read link:doc/for_developers.html[doc/for_developers.html] for an overview of how to write your own programs which build on Swarm-NG functionality.


Buliding the Swarm-NG Libraries and Demo Programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.What hardware & software is required to build & run Swarm-NG?
* Hardware
** Linux-based workstation (certain kernel versions?)
** nVidia GPU with G80 or G92 chip (we intend to support Fermi in near future)
* Software
** CUDA v2.3 and/or v3.0 
** g++ v?.? (later versions likely to work as well)
** Boost++ v?.? (later versions likely to work as well)
** Anything else?

.How do I build Swarm-NG for Linux? [Mario]
Run `make`.  Users can specify the paths to compiler, libaries, and include files, as well as compiler flags in the file ''Makefile.user''.  

.How do I build Swarm-NG for Windows? [Young In]
While we do not intend to support Swarm-NG for Windows, we can pass along some of the following tips.

* Supported host compiler on Windows platform is The Microsoft Visual Studio compiler, cl. 
* Supported build environments for nvcc on Windows are 
** Windows DOS shell
** Windows CygWin shells, use nvcc's drive prefix options.
** Windows MinGW shells, use nvcc's drive prefix options.

Although a variety of POSIX style shells is supported on Windows, nvcc will still
assume the Microsoft Visual Studio compiler for host compilation. Use of gcc is not
supported on Windows. You can find more information from `/usr/local/cuda/doc/nvcc_3.0.pdf`


.How do I test Swarm-NG? [Mario]
`make test`

Demo Programs [authors of each demo code]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.Commands to run demo codes
First, build the tutorials and demos with `make`.  Then, you can run the following demos from the tutorial.

- swarm_tutorial_cpu:  run `bin/swarm_tutorial_cpu`
- swarm_tutorial_gpu:  run `bin/swarm_tutorial_gpu`
- swarm_tutorial_compare:  run `bin/swarm_tutorial_compare`
- swarm_tutorial_montecarlo:  change into the run sub-directory and run `../bin/swarm_tutorial_montecarlo integrator-tutorial_montecarlo.cfg`
The behavior of this tutorial can be easily adjusted by changing the configuration file.  See the section below on xref:<DemoOptions>[Demo Options] for a list of required/valid parameters.

[[DemoOptions]]
.Demo Program Features/Options [authors of each demo program]

- Parameters for swarm_tutorial_montecarlo include:
** integrator: Specifies which integration algorithm will be used
** runon: Must be "cpu" or "gpu".  In v0.1, this must match the type of integrator specified by "integrator"
** dT: Specifies the duration of the integration
** h: Specifies the timestep of fixed timestep algorithms (e.g., hermite, rk4)
** Complete list

- Parameters for swarm include:
** integrator: Specifies which integration algorithm will be used
** runon: Must be "cpu" or "gpu".  In v0.1, this must match the type of integrator specified by "integrator"
** dT: Specifies the duration of the integration
** h: Specifies the timestep of fixed timestep algorithms (e.g., hermite, rk4)
** Complete list

.Known Limitations and/or Bugs [Young In]

- Limited to N-body systems with 3-10 bodies
- Must be many more!

Organization of Swarm-NG
~~~~~~~~~~~~~~~~~~~~~~~~

.Organization of files [Young In]
Swarm-NG has the following listing in the *swarm* directory. 
------------------------------
  <DIR>     docs              <1>
  <DIR>     integrators       <2>
  <DIR>     run               <3>
  <DIR>     scripts           <4>
  <DIR>     src               <5>
  <DIR>     test              <6>
            .gitignore
            Doxyfile
            Makefile
            README
            swarm.kdevelop
            swarm.kdevelop.filelist
------------------------------

<1> <<anchor-d1,docs>> contains documentations (tutorial, ...)
<2> <<anchor-d2,integrators>> contains integrators (verlet, hermite, ...)
<3> <<anchor-d3,run>> contains ...
<4> <<anchor-d4,scripts>> contains ...
<5> <<anchor-d5,src>> contains ...
<6> <<anchor-d6,test>> contains ...


[[anchor-d1]]
*docs*
------------------------------
            build_system.txt  <1>
            eventlog.txt      <2>
            README.swarm_adap <3>
            snapshotting.txt  <4>
------------------------------

<1> build system
<2> eventlog
<3> README
<4> snapshotting


[[anchor-d2]]
*integrators*
------------------------------
  <DIR>     euler                                      <1>
                 euler.cpp
                 euler.cu
                 euler.h
                 Makefile.mk

  <DIR>     hermite_adap_gpu                           <2>
                 hermite_adap_gpu.cpp
                 hermite_adap_gpu.cu
                 hermite_adap_gpu.h
                 hermite_adap_gpu_integrator_body.cu
                 Makefile.mk

  <DIR>     hermite_cpu                                <3>
                 hermite_cpu.cpp
                 hermite_cpu.h
                 Makefile.mk
                 ThreeVector.hpp

  <DIR>     hermite_gpu                                <4>
                 hermite_gpu.cpp
                 hermite_gpu.cu
                 hermite_gpu.h
                 hermite_gpu_integrator_body.cu
                 Makefile.mk

  <DIR>     mvs                                        <5>
                 Makefile.mk
                 mvs.cpp
                 mvs.cu
                 mvs.h

  <DIR>     rk4                                        <6>
                 Makefile.mk
                 rk4.cpp
                 rk4.cu
                 rk4.h

  <DIR>     verlet                                     <7>
                 Makefile.mk
                 verlet.cpp
                 verlet.cu
                 verlet.h

  <DIR>     verlet_cpu                                 <8>
                 Makefile.mk
                 verlet_cpu.cpp
                 verlet_cpu.h

------------------------------

<1> euler gpu
<2> hermite adaptive gpu
<3> hermite cpu
<4> hermite gpu
<5> mvs gpu
<6> rk4 gpu
<7> verlet gpu
<8> verlet cpu


[[anchor-d3]]
*run*
------------------------------
           HD_82943.rvin
           integrator-acboley.cfg
           integrator-adap.cfg
           integrator-mjuric.cfg
           integrator-tutorial_montecarlo.cfg
           integrator-tutorial_montecarlo_rv.cfg
           integrator.cfg
           README
------------------------------

[[anchor-d4]]
*scripts*
------------------------------
  <DIR>    ic
              peaShooterGen.py             <1>
              swarm_adap.py                <2>

           combine_cu_files.sh             <3>
           easyGen.py                      <4>
           generate_app_makefiles.sh       <5>
           peaShooterGen.py                <6>
           run_tests.sh                    <7>
           silence_nvcc_warnings.sh        <8>
           tarball_head.sh                 <9>
           throwBinaryGen.py               <10>
------------------------------
<1> ...
<2> ...
<3> ...
<4> ...
<5> ...
<6> ...
<7> ...
<8> ...
<9> ...
<10> ...


[[anchor-d5]]
*src*
------------------------------
  <DIR>    astro                                 <1>
                 BinaryStream.cpp
                 binarystream.h
                 constants.h
                 macros.h
                 MemoryMap.cpp
                 memorymap.h
                 types.h
                 Util.cpp
                 util.h

  <DIR>    cux                                   <2>
                 cuda_rng.h
                 cux.cpp
                 cux.h
                 cux_lowlevel.h

  <DIR>    gpulog                                <3>
  <DIR>          bits
                     gpulog_align.h
                     gpulog_constants.h
                     gpulog_debug.h
                     gpulog_ilogstream.h
                     gpulog_log.h
                     gpulog_logrecord.h
                     gpulog_macro_cleanup.h
                     gpulog_msg_layout.h
                     gpulog_printf.h
                     gpulog_ttraits.h
                     gpulog_types.h
                     gpulog_write.h

                 gpulog.h
                 lprintf.h

  <DIR>    scatter                               <4>
                 peaShooter.cpp

  <DIR>    swarm_adap                            <5>
                 COPYING.swarm_adap
                 swarm_adap.cpp
                 swarm_adap.h

           stopwatch.h                           <6>
           swarm.cpp                             <7>
           swarm.h                               <8>
           swarm_test_energy.cpp                 <9>
           swarm_tutorial_compare.cpp            <10>
           swarm_tutorial_cpu.cpp                <11>
           swarm_tutorial_gpu.cpp                <12>
           swarm_tutorial_montecarlo.cpp         <13>
           swarm_tutorial_montecarlo_rv.cpp      <14>
           swarmio.cpp                           <15>
           swarmio.h                             <16>
           swarmlib.cpp                          <17>
           swarmlib.cu                           <18>
           swarmlog.cpp                          <19>
           swarmlog.h                            <20>
           swarmquery.cpp                        <21>
           ThreeVector.hpp                       <22>
------------------------------

[[anchor-d6]]
*test*
------------------------------
  <DIR>    swarm
                filter
                integrator.cfg
                output.ref
                swarm-run.create
                swarm-run.test
------------------------------


.Overview of file contents [Young In]
Add text/lists here

.Where can I find additional documentation?
There is additional documentation for various componets of Swarm-NG in the docs subdirectory.
As the library interface solidifies, we anticipate providing additional documention in future releases and/or at the http://www.astro.ufl.edu/~eford/code/swarm/[Swarm-NG Website]

.What should I do if I think I found a compile-time bug?
- Please first check that you can successfully run `make test`.  If this fails, then please double check that you have the required drivers, libraries installed and paths set correctly.  You may need to create a Makefile.user to specify the path to your compilers, libraries, include files.  You may also need to specify the path to your cuda libraries via LD_LIBRARY_PATH.
- Please identify the minimal ammount of code that produces the compiler error.
- Contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list']

.What should I do if I think I found a run-time bug?
Unfortunately, GPUs can sometime get stuck in a problematic state, particularly during code development and when not running in console mode.  Before reporting a runtime bug, please test the code immediately after reboot your system into console mode.  
If the problem persists, then please check the archives of the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'] at the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG] to see if the bug has already been disussed.
If not, then please contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].

.How can I learn more about and/or contribute to the the Swarm-NG project?
- Visit the http://www.astro.ufl.edu/~eford/code/swarm/[Swarm-NG Website].
- Complete the http://astro.ufl.edu/~eford/code/swarm/[Swarm-NG User Questionnaire]
- Download, install and user Swarm-NG code (via link above)
- Read the archives of the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'] at the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG].
- Join the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG] to sign up for announcements via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].
- Contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].

.Licensing Information
Swarm-NG v0.1 may be licensed under the GPL v3, avaliable link:doc/gpl-3.0.txt[here].

.Acknowledging use of Swarm-NG
When writing papers which made use of Swarm-NG, we strongly encourage you to cite the technical paper descibing Swarm-NG.  We will provide a citation as soon as a preprint of this paper is avaliable.

.Acknowledgements
Contributors to the Swarm-NG codebase include: Eric Ford, Ameya Gancchha, Jianwei Gao, Mario Juric and Young In Yeo.  
We have adapted code from the GNU Science Library and John Chamber's Mercury code, which builds on a long history of n-body integrators, including Hal Levison and Martin Duncan's Swift code.
We greatfully acknowledge financial support from the University of Florida's Research Opportunity Seed Fund and the NASA Applied Information Systems Research Program.

Older Text (to be incorporated into this or other Documentation Files)
----------------------------------------------------------------------
== General layout ==

- swarm.cpp
        Main program file (main() is here)
- swarm.h
        Declarations of ensemble and integrator interfaces
- swarmlib.cpp
        Library of supporting host-only functions for the integrator, ensembles, etc.
- swarmlib.cu
        Library of supporting device functions for the integrators

- swarm.cu
        Autogenerated file (don't edit it by hand!) consisting of include
        statements including swarmlib.cu and all integ.*.cu files from the
        current directory.

== Adding integrators ==

To add a new integrator, follow the following simple steps:

- create a subdirectory integrators/<integrator_name>
- place any source files your integrator requires to that subdirectory
- create a Makefile.mk file there, based on example Makefile.mk files found with
  other integrators (e.g., look at integrators/euler/Makefile.mk)

Your integrator must:

- derive from 'integrator' class, and override one its virtual integrate()
  methods.
- provide an 'extern "C" integrator *create_XXXXX(const config &cfg)'
  function (where XXXX is the name of the integrator) that will return
  an instance of the integrator when called. Example:

        // factory
        extern "C" integrator *create_gpu_euler(const config &cfg)
        {
                return new gpu_euler_integrator(cfg);
        }

== Existing integrators: description and configuration parameters ==

* Hermite
        - h (real number)
                Timestep of the integrator
        - precision (integer)
                1 == double precision, 3 == mixed precision
