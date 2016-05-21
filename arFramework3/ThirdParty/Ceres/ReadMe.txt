This folder contains the compiler for the Matlab Interface of Ceres Nonlinear Solver. 

ceresd2d.cpp contains the mex interface for Ceres
compileCeres.m contains the compiling Script for ceresd2d.cpp

ceres-solver/ contains Ceres 1.11 as available at http://ceres-solver.org
eigen3/       contains Eigen 3.2.8 as availbe at http://eigen.tuxfamily.org/index.php?title=Main_Page


License 
    for Ceres can be found in ceres-solver/LICENSE
    for Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Licensing_FAQ.
    Eigen is under MPL2 license.


////////////////////////// CHANGES IN CERES-SOLVER FOLDER ///////////////////////

To enable compiling of the Ceres cc files without prior building on indivividual 
systems the Ceres config file config.h had to be manually created and inserted in
ceres-solver/include/ceres/internal/

When updating Ceres to a newer version be aware that this process has to be done 
manually and possibly further compiling options have to be added, depending on 
the operating system.



