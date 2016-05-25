This folder contains the compiler for the Matlab Interface of Ceres Nonlinear Solver. 

Please be aware that the compiling of Ceres requires a compiler with support for 
variadic templates:
Visual Studio C++   supported from version 2013
GCC                 supported from version 4.4
Clang C++           supported from version 2.9


ceresd2d.cpp    contains the mex interface for Ceres
ceresd2d.h      contains the relevant includes
compileCeres.m  contains the compiling Script for ceresd2d.cpp

ceres-solver/ contains Ceres 1.11   as available at http://ceres-solver.org
eigen3/       contains Eigen 3.2.8  as availbe at http://eigen.tuxfamily.org/index.php?title=Main_Page


License 
    for Ceres can be found in ceres-solver/LICENSE.
    Ceres is under New BSD license
    for Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Licensing_FAQ.
    Eigen is under MPL2 license.


////////////////////////// CHANGES IN CERES-SOLVER FOLDER ///////////////////////

To enable compiling of the Ceres cc files without prior building on indivividual 
systems the Ceres config file config.h had to be manually created and inserted in
ceres-solver/include/ceres/internal/

When updating Ceres to a newer version be aware that this process has to be done 
manually and possibly further compiling options have to be added, depending on 
the operating system.



//////////////////////// GCC VERSION ERROR WHEN RUNNING CERES ///////////////////

When compiling Ceres with a GCC compiler, whose version is above 4.7.x (the currently max
version supported by MATLAB), you might run into an error when running Ceres which looks like this:

.../libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by …)

This is due to Matlab preferably using the c++ standard library it has shipped with its distribution.
This libstdc++, located in MATLABROOT/sys/os/ is outdated and not the one your compiler used when 
compiling Ceres.

A fix for this issue is to symbolicly link MATLABs shipped version to the newer one your compiler used.
This can be done (on linux) via the commands:

// Backup old version
    mv MATLABROOT/sys/os/libstdc++.so.6 MATLABROOT/sys/os/libstdc++.so.6.old

// Create symbolic link to local libstdc++ folder
    ln -s /path/to/local/libstdc++.so.6 MATLABROOT/sys/os/libstdc++.so.6

If you don't know where your local libstdc++.so.6 is located just run
    locate libstdc++.so.6


////////////////////////////////////////////////////////////////////////////////////////////////////////
An alternate way of fixing this problem is to explicitly tell MATLAB to use the system libstd library.
However this has to be done before running matlab with the following command and is required for every run:
env LD_PRELOAD=/path/to/local/libstdc++.so.6 matlab -desktop





