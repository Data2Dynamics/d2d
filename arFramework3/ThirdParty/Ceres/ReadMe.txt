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



