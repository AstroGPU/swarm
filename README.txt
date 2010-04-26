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
* Read link:doc/build_system.html[doc/build_system.html] for more info on how to build the library and example programs.  
* Do a few quick tests with 'make test'.
* See how much of a speedup you get with 'make benchmark'.  (Warning, this requires performing the integrations on both the CPU and GPU, so it can take a long time.  Please be patient.)

* To learn how to access the Swarm-NG library from you own code, work through the series of introductory tutorials in src/tutorials (in the following order).  
. 'src/tutorials/swarm_tutorial_cpu.cpp': uses swarm library with CPU)
. 'src/tutorials/swarm_tutorial_gpu.cpp': uses swarm library with GPU
. 'src/tutorials/swarm_tutorial_compare.cpp': compares results from GPU & CPU
. 'src/tutorials/swarm_tutorial_benchmark.cpp': compares performance using GPU & CPU
. 'src/tutorials/swarm_tutorial_montecarlo.cpp': performs many n-body integrations in Monte Carlo fashion in user-specified region of parameter space.  This program uses a configuration file.  The format is desribed in link:docs/swarm_tutorial_montecarlo.html[docs/swarm_tutorial_montecarlo.html].
. 'src/swarm.cpp': this reads in an ensemble of initial conditions and writes output to binary files.  The output can be easily accessed via the swarmquery program.

* To see how an end user can use 'swarm' as is:
** Take a look at the test 'test/swarm/swarm-run.test'.  If you run 'make test', the output will be stored in 'test-outputs/swarm/swarm-run.test/output.log'
** Read link:doc/configuration_file.html[doc/configuration_file.html] for a desription of the configuration file format and valid parameters for the 'bin/swarm' program.  (Some of the demos use a configuration file with the same format, but they may not accept all the parameters.)
** To see how end users can access data written by the swarm logging system, use 'swarmquery output_file'.  For info on options, run 'swarmquery --help' 

* To learn how developers an use the swarm library:
** Read link:doc/for_developers.html[doc/for_developers.html] for an overview of how to write your own programs which build on Swarm-NG functionality.
** Take a look at 'doc/swarm_scatter_demo.html' and 'src/scatter_demo/swarm_scatter_demo.cpp' for an example 
** Read link:doc/eventlog.txt[doc/eventlog.txt] and link:snapshotting.txt[do/snapshotting.txt] for info on using the swarm output system.

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
** Asciidoc & Doxygen (if you want to generate html documentation)
** Anything else?

.How do I build Swarm-NG for Linux? 
Run `make`.  Users can specify the paths to compiler, libaries, and include files, as well as compiler flags in the file ''Makefile.user''.  

.How do I build Swarm-NG for Windows? 
While we do not intend to support Swarm-NG for Windows, we can pass along some of the following tips.

* We're told that following build environments for nvcc work on Windows using the The Microsoft Visual Studio compiler, cl, as the host compiler:
** Windows DOS shell
** Windows CygWin shells, use nvcc's drive prefix options.
** Windows MinGW shells, use nvcc's drive prefix options.

* Although a variety of POSIX style shells is supported on Windows, nvcc will still assume the Microsoft Visual Studio compiler for host compilation. Use of gcc is not supported on Windows. 

* You can find more information from `/usr/local/cuda/doc/nvcc_3.0.pdf`

.How do I test Swarm-NG? 
`make test`


Demo Programs 
~~~~~~~~~~~~~

.Commands to run demo codes
First, build the tutorials and demos with `make`.  Then, you can run the following demos from the tutorial.

- swarm_tutorial_cpu:  run `bin/swarm_tutorial_cpu`
- swarm_tutorial_gpu:  run `bin/swarm_tutorial_gpu`
- swarm_tutorial_compare:  run `bin/swarm_tutorial_compare`
- swarm_tutorial_benchmark:  run `bin/swarm_tutorial_benchmark`.  Several options are avaliable.  run with the command line option --help for more info.
- swarm_tutorial_montecarlo:  change into the run sub-directory and run `../bin/swarm_tutorial_montecarlo integrator-tutorial_montecarlo.cfg`
The behavior of this tutorial can be easily adjusted by changing the configuration file.  See docs/swarm_tutorial_montecarlo.html[docs/swarm_tutorial_montecarlo.html] for details.
- swarm_scatter_demo:  This requires multiple input files.  See link:docs/swarm_scatter_demo.html[docs/swarm_scatter_demo.html] for details.
- 
the section below on xref:<DemoOptions>[Demo Options] for a list of required/valid parameters.

.Known Limitations and/or Bugs [Young In]

- swarm is limited to N-body systems with 3-10 bodies.
- swarm only integrates forward in time.
- While there are euler and rk4 integrators, they are in a preliminary form and should not be used.  
- E.g., rk4 requires exactly 3 bodies and does not include a cpu version

Organization of Swarm-NG
~~~~~~~~~~~~~~~~~~~~~~~~

.Organization of documentation
In addition to this README.txt in the root directory, there are several text files (ending in .txt or .man) in the *docs* subdirectory. +
If you have asciidoc installed, then you can generate pretty html (or other formats) with 'make doc-asciidoc'. +
While documenting source code is an ongoing projet, there is significant inline documentation.  If you have doxygen installed, then you can generate pretty html by running 'make doc-doxygen'. It places output in the *reference* subdirectory.+

.Organization of source code
The main source code is in the *src* directory. +
Source code for the tutorials is in the *src/tutorials* directory. + 
Source code for the swarm_scatter_demo is in the *src/scatter_demo* directory. + 
Source code for the various integrators is in subdirectories of *src/integrators*. + 
Source code for the GPU logging system is in the *src/gpulog* directory. +
Various utilities are in *src/cux* and *src/astro*. +

.Organization of executables
All swarm executables and the swarm libary (libswarm.so) are placed in the *bin* directory. +

.Other directories
The *scripts* directory has scripts needed for the build system, as well as for generating initial conditions for some demo programs.  +
The *run* directory contains configuration files and is where swarm is intended to be run from. 

For more information on the organization of directories and files, see link:docs/organization.txt[docs/organization.txt] and/or the link:reference/html/files.html[doxygen files page].

Digging Deeper
~~~~~~~~~~~~~~
.Where can I find additional documentation?
There is additional documentation for various componets of Swarm-NG in the docs subdirectory.  Most are written in asciidoc format.  To generate html versions, run 'make doc-asciidoc'.+
After running 'make doc-doxygen', you can also browse the link:reference/html/index.html[doxygen reference]. +

As the library interface solidifies, we anticipate providing additional documention in future releases and/or at the http://www.astro.ufl.edu/~eford/code/swarm/[Swarm-NG Website]

.What should I do if I think I found a compile-time bug?
- Please first check that you can successfully run `make test`.  If this fails, then please double check that you have the required drivers, libraries installed and paths set correctly.  You may need to create a Makefile.user to specify the path to your compilers, libraries, include files.  You may also need to specify the path to your cuda libraries via LD_LIBRARY_PATH.
- If the error arises while compiling the included tests, then run 'make feedback', whih will generate a file named 'feedback.tgz' which we may request that you send us.
- Otherwise, please identify the minimal ammount of code that produces the compiler error.
- Contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list']

.What should I do if I think I found a run-time bug?
Unfortunately, GPUs can sometime get stuck in a problematic state, particularly during code development and when not running in console mode.  Before reporting a runtime bug, _please test the code immediately after reboot your system into console mode_.  
If the problem persists, then please check the archives of the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'] at the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG] to see if the bug has already been disussed.
If not, then please contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].

.How can I learn more about and/or contribute to the the Swarm-NG project?
- Visit the http://www.astro.ufl.edu/~eford/code/swarm/[Swarm-NG Website].
- Complete the http://astro.ufl.edu/~eford/code/swarm/[Swarm-NG User Questionnaire]
- Download, install and user Swarm-NG code (via link above)
- Read the archives of the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'] at the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG].
- Join the http://groups.google.com/group/swarm-ng[Google Group for Swarm-NG] to sign up for announcements via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].
- Contact Swarm-NG developers via the mailto:swarm-ng@googlegroups.com['Swarm-NG mailing list'].

Licensing, Citing, Acknowledgements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.Licensing Information

Swarm-NG v0.1 may be licensed under the GPL v3, avaliable
link:doc/gpl-3.0.txt[here].

.Acknowledging use of Swarm-NG

When writing papers which made use of Swarm-NG, we strongly encourage
you to cite the technical paper descibing Swarm-NG.  We will provide a
citation as soon as a preprint of this paper is avaliable.

.Acknowledgements Swarm-NG is the product of a collaboration between
astrophysicists (UF Astronomy department and Harvard-Smithsonian
Center for Astrophysics) and computer scientists at UF's
http://www.cise.ufl.edu/research/SurfLab/[SurfLab].  Contributors to
the Swarm-NG codebase include: Eric Ford, Ameya Gancchha, Jianwei Gao,
Mario Juric and Young In Yeo.  +

We have adapted code from the the nVidia CUDA SDK, the GNU Science
Library, Sverre Aarseth's hermit2 code, and John Chamber's Mercury
code, which builds on a long history of n-body integrators, including
Hal Levison and Martin Duncan's Swift code. +

We greatfully acknowledge financial support from the University of
Florida's Research Opportunity Seed Fund, the National Science
Foundation (grant CCF-0728797) and especially the NASA Applied
Information Systems Research Program (grant NNX09AM41G).  We also
acknowledge support from the UF High Performance Computing Center and
nVidia Corporation.



