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

The software was awarded twice as **best performer** in the *Dialogue for Reverse Engineering Assessments and Methods* ([DREAM](http://www.the-dream-project.org/)). 2011 in the *Estimation of Model Parameters Challenge* ([DREAM6](http://www.the-dream-project.org/challenges/dream6-estimation-model-parameters-challenge)) and 2012 in the *Network Topology and Parameter Inference Challenge* ([DREAM7](http://www.the-dream-project.org/challenges/network-topology-and-parameter-inference-challenge)). Read more about this in: Steiert B., et al. [Experimental Design for Parameter Estimation of Gene Regulatory Networks](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0040052). PLoS ONE 7(7), e40052, 2012

## Further topics

Please read the [Wiki](https://github.com/Data2Dynamics/d2d/wiki/Home/) if you are interested in details.


## Copyright notice
Copyright © 2016 D2D development team. All rights reserved. [Copyright text and third party software license information.](https://github.com/Data2Dynamics/d2d/wiki/Copyright)
