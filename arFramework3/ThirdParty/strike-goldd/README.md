# About STRIKE-GOLDD

STRIKE-GOLDD is a MATLAB toolbox that analyses nonlinear models of ordinary differential equations. It performs a simultaneous assessment of:
- state **observability**,
- parameter structural local **identifiability**,  
- unknown **input observability** (supported since v2.1.1). 

The analysis is performed symbolically, and yields results that are valid for all values of the variables, except for a set of measure zero. More information is available at the [STRIKE-GOLDD website](https://sites.google.com/site/strikegolddtoolbox/home).

STRIKE-GOLDD was created by [Alejandro F. Villaverde](https://sites.google.com/site/alexfvillaverde/), <afvillaverde@uvigo.gal>. 

## Installation and requirements

STRIKE-GOLDD requires a MATLAB installation with the Symbolic toolbox. Additionally, if one wishes to use optimization-based decomposition, the MATLAB version of the [MEIGO](http://nautilus.iim.csic.es/~gingproc/meigo.html) toolbox is required.

To use STRIKE-GOLDD you just need to:
1. download the code
2. open a MATLAB session
3. run `install` from the STRIKE-GOLDD directory
4. define the problem by editing `options.m`
5. run `STRIKE-GOLDD`

More information can be found in the [STRIKE-GOLDD manual](STRIKE-GOLDD/doc/STRIKE-GOLDD_manual.pdf)

## Publications

Publication of the methodology (first version of STRIKE-GOLDD):

[Villaverde AF, Barreiro A, Papachristodoulou A (2016). Structural Identifiability of Dynamic Systems Biology Models. *PLoS Computational Biology* 12(10):e1005153, doi:10.1371/journal.pcbi.1005153](http:dx.doi.org/doi:10.1371/journal.pcbi.1005153)

Extension for time-varying inputs (STRIKE-GOLDD 2):

[Villaverde AF, Evans ND, Chappell MJ, Banga JR (2019). Input-dependent structural identifiability of nonlinear systems. *IEEE Control Systems Letters* 3(2):1–6, doi:10.1109/LCSYS.2018.2868608](http://dx.doi.org/doi:10.1109/LCSYS.2018.2868608)

Extension for unknown inputs; FISPO analysis (STRIKE-GOLDD 2.1):

[Villaverde AF, Tsiantis N, Banga JR (2019). Full observability and estimation of unknown inputs, states, and parameters of nonlinear biological models. *Journal of the Royal Society Interface* 16(156), doi:10.1098/rsif.2019.0043](http://dx.doi.org/doi:10.1098/rsif.2019.0043)

Extension for finding Lie symmetries (STRIKE-GOLDD 2.1.6):

[Massonis G & Villaverde AF (2020). Finding and breaking Lie symmetries: implications for structural identifiability and observability in biological modelling. *Symmetry* 12(3), 469, doi:10.3390/sym12030469](https://doi.org/10.3390/sym12030469)

Extension for multi-experiment analysis and implementation of the ORC-DF algorithm (STRIKE-GOLDD 2.2):

[Martínez N & Villaverde AF (2020). Nonlinear observability algorithms with known and unknown inputs: analysis and implementation. *Mathematics* 8(11), 1876, doi:10.3390/math8111876](https://doi.org/10.3390/math8111876)

Extension for automatic reparameterization, AutoRepar (STRIKE-GOLDD 3.0):

[Massonis G, Banga JR, Villaverde AF (2020). Repairing dynamic models: a method to obtain identifiable and observable reparameterizations with mechanistic insights. *arXiv*, arXiv:2012.09826 [eess.SY]](https://arxiv.org/abs/2012.09826)


## Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
    
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
