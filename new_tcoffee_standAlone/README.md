BioTools++
==========

1. Requirements
---------------

The following libraries need to be installed:
- CMAKE
- OpenMP
- Boost (>1.46)

In case the Boost library is not installed in the system path please set the
environment variable BOOST_ROOT to the correct directory. 

For example:
export BOOST_ROOT=/usr/soft/boost 


2. Installation
---------------


You can either compile a single program you choose (a) or you compile and install all of them at once (b)

(a) Compile a single program of this suit:

1. cd build
2. cmake ..
3. make <program name>

(b) Compile and install all of them:

1. cd build
2. cmake -DCMAKE_INSTALL_PREFIX=<path of installation directory> ..
3. make install

If you use cmake as in (a) but without the "program name" all programs will be compiled but not copied to another directory.

3. The manual
-------------
The manuals for the different programs can be found in the doc/manual directory.


4. Documentation
----------------

The function documentation with doxygen can be generated after the cmake configuration using the following
command:

make doc


5. Problems
-----------

This is still a beta version. In case you encounter any problem please contact me.
