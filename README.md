# Data2Dynamics Software

**Contact:** Andreas Raue - <andreas.raue@fdm.uni-freiburg.de>
(for support issues, please use the issues and forum, thanks!)

**Cite:** 

* [Data2Dynamics: a modeling environment tailored to parameter estimation in dynamical systems.](https://academic.oup.com/bioinformatics/article/31/21/3558/195191) Raue A., et al. Bioinformatics, 31(21), 3558-3560, 2015.

* [Lessons Learned from Quantitative Dynamical Modeling in Systems Biology.](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0074335) Raue A., et al. PLOS ONE, 8(9), e74335, 2013.

**Video:** See how easy model calibration can be using D2D Software in this [YouTube Video](http://www.youtube.com/watch?v=_aAtSo_xe7I). The video is recorded in real time on a dual core laptop with the [Bachmann et al. 2011](https://github.com/Data2Dynamics/d2d/wiki/Bachmann_MSB2011) example using the first data set.

## Features

Major purpose: Establishing ODE models based on experimental data. The software is designed for biochemical reaction networks, but is not limited to this. 

Key feature: Reliable and efficient parameter estimation techniques and statistical assessment of parameter-, measurement- and prediction uncertainties.

Some special features:

* The framework can deal with xperimental error bars but also allows fitting of error parameters (error models).
 
* Model inputs can be implemented as parameterized functions or cubic splines and can be estimated together with the model dynamics ([read more](https://github.com/Data2Dynamics/d2d/wiki/Input-estimation)).

* For model calibration, i.e. parameter estimation, both stochastic and deterministic numerical optimization algorithms can be used.

* For uncertainty analysis of model parameter and predictions, the profile likelihood approach and Markov chain Monte Carlo sampling approaches are available ([read more](https://github.com/Data2Dynamics/d2d/wiki/Uncertainty-analysis)).

* Efficient numerical calculation of the dynamics and derivatives in a parallelized manner ([read more](https://github.com/Data2Dynamics/d2d/wiki/Parallelization)). 

* L1 regularization of parameter fold-changes can be used ([read more](https://github.com/Data2Dynamics/d2d/wiki/L1-regularization)).

* L2 regularization and incorporation of prior knowledge 

* Identification of informative experimental designs 

* The software is open source and free for non-commercial usage and supports parallelization

A more comprehensive description of features is available [features page](https://github.com/Data2Dynamics/d2d/wiki/Features) in the [Wiki](https://github.com/Data2Dynamics/d2d/wiki/Home).

## Awards

The software was awarded twice as **best performer** in the *Dialogue for Reverse Engineering Assessments and Methods* ([DREAM](http://www.the-dream-project.org/)). 2011 in the *Estimation of Model Parameters Challenge* ([DREAM6](http://www.the-dream-project.org/challenges/dream6-estimation-model-parameters-challenge)) and 2012 in the *Network Topology and Parameter Inference Challenge* ([DREAM7](http://www.the-dream-project.org/challenges/network-topology-and-parameter-inference-challenge)). Read more about this in: Steiert B., et al. [Experimental Design for Parameter Estimation of Gene Regulatory Networks](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0040052). PLoS ONE 7(7), e40052, 2012

## Further topics

Please read the [Wiki](https://github.com/Data2Dynamics/d2d/wiki/Home/) if you are interested in details.


## Copyright notice
Copyright © 2016 D2D development team. All rights reserved. [Copyright text and third party software license information.](https://github.com/Data2Dynamics/d2d/wiki/Copyright)
