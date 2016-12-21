Original folder name: E:\clemens\systemBiologie_cvs\Projekte\B3-InsulinIgf1\2014_Ole\2_MitElisa 
Date: 21-Dec-2016 

Model states and right hand side of the ODEs:
Insulin		BoundUnspec*koff_unspec + IR1*kd1 - Ins*kon_unspec - Ins*Rec1*ka1 + IR2*kd1*kd2fold - Ins*Rec2*ka1*ka2fold
Receptors_low		IR1*kd1 + IR1in*kout_frag - Ins*Rec1*ka1
Receptors_high		IR2in*kout_frag + IR2*kd1*kd2fold - Ins*Rec2*ka1*ka2fold
IR_low		IR1in*kout - IR1*kin - IR1*kd1 + Ins*Rec1*ka1
IR_high		IR2in*kout2 - IR2*kin2 - IR2*kd1*kd2fold + Ins*Rec2*ka1*ka2fold
IR_low_internalized		IR1*kin - IR1in*kout - IR1in*kout_frag
IR_high_internalized		IR2*kin2 - IR2in*kout2 - IR2in*kout_frag
Uptake1		Ins*Rec1*ka1 - IR1*kd1
Uptake2		Ins*Rec2*ka1*ka2fold - IR2*kd1*kd2fold
InsulinFragments		IR1in*kout_frag + IR2in*kout_frag
BoundUnspec		Ins*kon_unspec - BoundUnspec*koff_unspec

Observables and Errors:
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
IR1_obs:		fy=scale*(IR1 + offset + IR1in)		fystd=IR_obs_std
IR2_obs:		fy=scale*(IR2 + offset + IR2in)		fystd=IR_obs_std
IRsum_obs:		fy=scale*(0.605*(IR1+IR1in) + (1-0.605)*(IR2+IR2in) + offset)		fystd=IR_obs_std
Insulin_obs:		fy=offset_nExpID1 + scaleElisa_nExpID1*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID1)		fystd=std
Insulin_obs:		fy=offset_nExpID1 + scaleElisa_nExpID1*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID1)		fystd=std
Insulin_obs:		fy=offset_nExpID1 + scaleElisa_nExpID1*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID1)		fystd=std
Insulin_obs:		fy=offset_nExpID2 + scaleElisa_nExpID2*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID2)		fystd=std
Insulin_obs:		fy=offset_nExpID2 + scaleElisa_nExpID2*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID2)		fystd=std
Insulin_obs:		fy=offset_nExpID2 + scaleElisa_nExpID2*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID2)		fystd=std
Insulin_obs:		fy=offset_nExpID3 + scaleElisa_nExpID3*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID3)		fystd=std
Insulin_obs:		fy=offset_nExpID3 + scaleElisa_nExpID3*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID3)		fystd=std
Insulin_obs:		fy=offset_nExpID3 + scaleElisa_nExpID3*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID3)		fystd=std
Insulin_obs:		fy=offset_nExpID4 + scaleElisa_nExpID4*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID4)		fystd=std
Insulin_obs:		fy=offset_nExpID4 + scaleElisa_nExpID4*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID4)		fystd=std
Insulin_obs:		fy=offset_nExpID4 + scaleElisa_nExpID4*(Ins + fragments * InsulinFragments)/(1+(Ins + fragments * InsulinFragments)/km_nExpID4)		fystd=std


Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model

           name                 lb       value       ub          10^value        fitted   prior
#   1|  E| IR_obs_std         |       -5       -1.3         -1 | 1      +0.05 |       1 | uniform(-5,-1) 
#   2|   | fragments          |       +0         +1         +1 | 0         +1 |       1 | uniform(0,1) 
#   3|DI | ini_R1             |       -5       +1.8         +3 | 1        +59 |       1 | normal(1.80977,3^2) 
#   4|DI | ini_R2fold         |       -5       +1.2         +3 | 1        +16 |       1 | normal(1.19692,3^2) 
#   5|D  | ka1                |       -5       -2.4         +3 | 1    +0.0038 |       1 | normal(-4.00677,3^2) 
#   6|D  | ka2fold            |       -5      +0.69         +3 | 1         +5 |       1 | normal(-0.0186139,3^2) 
#   7|D  | kd1                |       -5      +0.96         +3 | 1       +9.1 |       1 | normal(-0.917793,3^2) 
#   8|D  | kd2fold            |       -5      +0.85         +3 | 1         +7 |       1 | normal(-0.105609,3^2) 
#   9|D  | kin                |       -5      -0.43         +3 | 1      +0.37 |       1 | uniform(-5,3) 
#  10|D  | kin2               |       -5      -0.27         +3 | 1      +0.54 |       1 | uniform(-5,3) 
     |   |                    |                                |              |         |      
#  11|   | km_nExpID1         |       +3         +8         +8 | 1     +1e+08 |       1 | uniform(3,8) 
#  12|   | km_nExpID2         |       +3         +8         +8 | 1     +1e+08 |       1 | uniform(3,8) 
#  13|   | km_nExpID3         |       +3         +8         +8 | 1     +1e+08 |       1 | uniform(3,8) 
#  14|   | km_nExpID4         |       +3         +8         +8 | 1     +1e+08 |       1 | uniform(3,8) 
#  15|D  | koff_unspec        |       -5         +1         +3 | 1        +10 |       1 | uniform(-5,3) 
#  16|D  | kon_unspec         |       -5       +1.3         +3 | 1        +20 |       1 | uniform(-5,3) 
#  17|D  | kout               |       -5       -1.3         +3 | 1     +0.045 |       1 | uniform(-5,3) 
#  18|D  | kout2              |       -5       -1.4         +3 | 1     +0.036 |       1 | uniform(-5,3) 
#  19|D  | kout_frag          |       -5         -2         +3 | 1     +0.011 |       1 | uniform(-5,3) 
#  20|   | offset             |       -5       +1.2         +3 | 1        +15 |       1 | uniform(-5,3) 
     |   |                    |                                |              |         |      
#  21|   | offset_nExpID1     |       -5       -1.6         +5 | 1     +0.025 |       1 | uniform(-5,5) 
#  22|   | offset_nExpID2     |       -5       -1.7         +5 | 1     +0.022 |       1 | uniform(-5,5) 
#  23|   | offset_nExpID3     |       -5       -2.1         +5 | 1    +0.0073 |       1 | uniform(-5,5) 
#  24|   | offset_nExpID4     |       -5      +0.19         +5 | 1       +1.6 |       1 | uniform(-5,5) 
#  25|   | scale              |       -5      -0.88         +3 | 1      +0.13 |       1 | uniform(-5,3) 
#  26|   | scaleElisa_nExpID1 |     +0.1      +0.54         +1 | 0      +0.54 |       1 | uniform(0.1,1) 
#  27|   | scaleElisa_nExpID2 |     +0.1      +0.57         +1 | 0      +0.57 |       1 | uniform(0.1,1) 
#  28|   | scaleElisa_nExpID3 |     +0.1       +0.1         +1 | 0       +0.1 |       1 | uniform(0.1,1) 
#  29|   | scaleElisa_nExpID4 |     +0.1         +1         +1 | 0         +1 |       1 | uniform(0.1,1) 
#  30|  E| std                |       -5      -0.58         +3 | 1      +0.26 |       1 | uniform(-5,3) 


Short documentation:

This is the Insulin-Receptor model as used in [1] to as part of a multiscale model.

The FACS data is the preprocessed result of a large amount of FACS. A comprehensive data processing pipeline
has been applied as described in [2]
The result of this data preprocessing is stored in Data/datfit.mat and translated to the D2D format in FacsData_unlog10.xlsx

In comparison to [2], a 2nd data type namely ELISA of Insulin depetion in the medium was available (four experiments). 
The estimates of six dynamic parameters obtained in [2] were used as weak priors.

Please not that five parameters are at the bounds. This could be resolved by exanding the bounds in the following way:
The bounds of the scale- and km parameters could be enlarged  via: 
ar.ub(11:14)=10;
ar.ub(26:29)=3;
Then the profiles looks better and LHS also looks fine.

However, a new global optimium emerges which is diffcult to be found. It's even not directly seen by LHS(1000). 
The new global optimum has fragments parameter close to zero and does not catch the observed dynamics of insulin binding and over-estimate the observation errors.
In fact, prior knowledge about the SD of the IR-binding data is available in Data/datfit.mat and could be utilized.


Literature:
[1] L.O. Schwen, A. Schenk, C. Kreutz J. Timmer, M.M. Bartolome Rodriguez, L. Kuepfer, T. Preusser
    Representative sinusoids for hepatic four-scale pharmacokinetics simulations, Plos One (2015) 10, e0133653 
[2] Kreutz, Statistical Approaches for Molecular and Systems Biology, PhD Thesis, University of Freiburg, 2011.
