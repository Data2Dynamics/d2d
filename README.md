# Data2Dynamics Software

**Contact:** Andreas Raue - <andreas.raue@fdm.uni-freiburg.de>

**Cite:** 

* [Data2Dynamics: a modeling environment tailored to parameter estimation in dynamical systems.](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btv405?ijkey=YPsnNzFC4CIzy5g&keytype=ref) Raue A., et al. Bioinformatics, 31(21), 3558-3560, 2015.

* [Lessons Learned from Quantitative Dynamical Modeling in Systems Biology.](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0074335) Raue A., et al. PLOS ONE, 8(9), e74335, 2013.

**Video:** See how easy model calibration can be using D2D Software in this [YouTube Video](http://www.youtube.com/watch?v=_aAtSo_xe7I). The video is recorded in real time on a dual core laptop with the [Bachmann et al. 2011](https://github.com/Data2Dynamics/d2d/wiki/Bachmann_MSB2011) example using the first data set.

## Features

The Data2Dynamics software package is a collection of numerical methods for quantitative dynamic modeling and a comprehensive model and data description language. The software facilitates the construction of dynamical models for biochemical reaction networks, but is not limited to this. Its key features are reliable and efficient model calibration and parameter estimation techniques using numerical optimization methods and the assessment of measurement and model uncertainties. To this end various approaches are provided.

* Measurement noise of experimental data can either be explicitly provided or simultaneously estimated with the model dynamics. For the latter approach a parameterized noise model is used.

* Model inputs can be implemented as parameterized functions or cubic splines and can be estimated together with the model dynamics ([read more](https://github.com/Data2Dynamics/d2d/wiki/Input estimation)).

* For model calibration, i.e. parameter estimation, both stochastic and deterministic numerical optimization algorithms can be used.

* For uncertainty analysis of model parameter and predictions, the profile likelihood approach and Markov chain Monte Carlo sampling approaches are available ([read more](https://github.com/Data2Dynamics/d2d/wiki/Uncertainty analysis)).

* For efficient numerical solution of the dynamics and derivate calculations parallelized numerical solvers are implemented ([read more](https://github.com/Data2Dynamics/d2d/wiki/Parallelization)). In addition some functions take advantage of PARFOR loops that are provided by the MATLAB Distributed Computing Toolbox ([read more](https://github.com/Data2Dynamics/d2d/wiki/Distributed Computing)).

* For identification of cell type-specific parameters, L1 regularization of parameter fold-changes can be used ([read more](https://github.com/Data2Dynamics/d2d/wiki/L1 regularization)).

* L2 regularization and incorporation of prior knowledge by penalized maximum likehood / maximum a posteriori estimation

* Identification of informative experimental designs

* The software is open source and free for non-commercial usage.

## Example Applications

[Brief summary of all available example models.](https://github.com/Data2Dynamics/d2d/wiki/Example models in D2D)

The software and the theoretical concepts implemented therein were developed based on Systems Biology applications. Some of the applications are provided as examples and blue prints for new applications:

* [IL-13 induced JAK2/STAT5 signaling model in MedB-1 cells](https://github.com/Data2Dynamics/d2d/wiki/Raia_CancerResearch2011) - Raia et al., Cancer Research 2010 

* [Epo receptor model in BaF3 cells](https://github.com/Data2Dynamics/d2d/wiki/Becker_Science2010) - Becker et al., Science 2010

* [Epo induced JAK2/STAT5 signaling model in CFU-E cells](https://github.com/Data2Dynamics/d2d/wiki/Bachmann_MSB2011) - Bachmann et al., Molecular Systems Biology 2011

* [Epo induced JAK2/STAT5 signaling model in CFU-E and H838 cells](https://github.com/Data2Dynamics/d2d/wiki/Merkle_PCB2016) - Merkle et al., PLOS Computational Biology 2016

* [STAT5 activation model in BaF3-EpoR cells](https://github.com/Data2Dynamics/d2d/wiki/Boehm_JProteomeRes2014) - Boehm et al., Journal of proteome research 2014

* [JAK/STAT signaling model](https://github.com/Data2Dynamics/d2d/wiki/Swameye_PNAS2003) - Swameye et al., PNAS 2003

* [Indigoidine synthesis in *E. coli*](https://github.com/Data2Dynamics/d2d/wiki/Beer_MolBioSyst2014) - Beer et al., Molecular BioSystems 2014

The software was awarded twice as **best performer** in the *Dialogue for Reverse Engineering Assessments and Methods* ([DREAM](http://www.the-dream-project.org/)). 2011 in the *Estimation of Model Parameters Challenge* ([DREAM6](http://www.the-dream-project.org/challenges/dream6-estimation-model-parameters-challenge)) and 2012 in the *Network Topology and Parameter Inference Challenge* ([DREAM7](http://www.the-dream-project.org/challenges/network-topology-and-parameter-inference-challenge)). Read more about this in: Steiert B., et al. [Experimental Design for Parameter Estimation of Gene Regulatory Networks](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0040052). PLoS ONE 7(7), e40052, 2012

## Further topics
* [Installation and system requirements](https://github.com/Data2Dynamics/d2d/wiki/Installation)
* [Setting up models](https://github.com/Data2Dynamics/d2d/wiki/Setting up models)
* [First steps](https://github.com/Data2Dynamics/d2d/wiki/First steps)
* [Advanced events and pre-equilibration](https://github.com/Data2Dynamics/d2d/wiki/Advanced events and pre-equilibration)
* [Computation of integration-based prediction bands](https://github.com/Data2Dynamics/d2d/wiki/Computation of integration-based prediction bands)

## Frequently asked questions
* [How is the architecture of the code and the most important commands?](https://github.com/Data2Dynamics/d2d/wiki/CodeArchitecture)
* [What are the most important fields of the global variable ar?](https://github.com/Data2Dynamics/d2d/wiki/What are the most important fields of the variable ar%3F)
* [What are the most important functions?](https://github.com/Data2Dynamics/d2d/wiki/What are the most important functions in d2d%3F)
* [Optimization algorithms available in the d2d-framework](https://github.com/Data2Dynamics/d2d/wiki/Optimization algorithms available in the d2d-framework)
* [Objective function, likelhood and chi-square in the d2d framework](https://github.com/Data2Dynamics/d2d/wiki/Objective function, likelhood and chi-square in the d2d framework)
* [How to set up priors?](https://github.com/Data2Dynamics/d2d/wiki/Priors in the d2d framework)
* [How to set up steady state constraints?](https://github.com/Data2Dynamics/d2d/wiki/How to set up steady state constraints) 
* [How do I restart the solver upon a step input?](https://github.com/Data2Dynamics/d2d/wiki/How do I reset the solver upon a step input)
* [How to deal with integrator tolerances?](https://github.com/Data2Dynamics/d2d/wiki/Integration tolerances and modification possibilities)
* [How to implement a bolus injection?](https://github.com/Data2Dynamics/d2d/wiki/How to implement a bolus injection) 
* [How to implement washing and an injection?](https://github.com/Data2Dynamics/d2d/wiki/How to implement washing and bolus injection) 
* [How to implement a moment ODE model?](https://github.com/Data2Dynamics/d2d/wiki/How to implement a moment ODE model) 
* [How to run PLE calculations on a Cluster?](https://github.com/Data2Dynamics/d2d/wiki/How to run PLE calculations on a Cluster)

## Copyright notice
Copyright © 2016 D2D development team. All rights reserved. [Copyright text and third party software license information.](https://github.com/Data2Dynamics/d2d/wiki/Copyright)
