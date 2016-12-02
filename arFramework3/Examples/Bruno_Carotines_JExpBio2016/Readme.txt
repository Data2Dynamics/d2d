The model has been published in:

Mark Bruno, Julian Koschmieder, Florian Wuest, Patrick Schaub, Mirjam Fehling-Kaschek, Jens Timmer, Peter Beyer, and Salim Al-Babili 
Enzymatic study on AtCCD4 and AtCCD7 and their potential to form acyclic regulatory metabolites
J. Exp. Bot. first published online October 6, 2016 doi:10.1093/jxb/erw356 

Date of creating this example: 02-Dec-2016 


Model states and right hand side of the ODEs:
beta-carotin		-bcar*kb1
cry		- bcry*kc1 - bcry*kc2
beta-10		bcar*kb1 - b10*kb2 + bcry*kc1
beta-io		b10*kb2 + bcar*kb1 + bcry*kc2
OH-beta-10		bcry*kc2 - kc4*ohb10 + k5*zea
OH-beta-io		bcry*kc1 + kc4*ohb10 + k5*zea
zea		-k5*zea


Observables and Errors:
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1
obcar:		fy=bcar		fystd=1
obcry:		fy=bcry		fystd=1
obio:		fy=bio		fystd=1
ob10:		fy=b10		fystd=1
oohb10:		fy=ohb10		fystd=1
ozea:		fy=zea		fystd=1


Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model

           name         lb       value       ub          10^value        fitted   prior
#   1|DI | init_b10   |       -5      +0.63         +3 | 1       +4.2 |       1 | uniform(-5,3) 
#   2|DI | init_bcar1 |       -5      +0.65         +3 | 1       +4.4 |       1 | uniform(-5,3) 
#   3|DI | init_bcar2 |       -5      +0.52         +3 | 1       +3.3 |       1 | uniform(-5,3) 
#   4|DI | init_bcry  |       -5      +0.72         +3 | 1       +5.2 |       1 | uniform(-5,3) 
#   5|DI | init_ohb10 |       -5      +0.24         +3 | 1       +1.7 |       1 | uniform(-5,3) 
#   6|DI | init_zea   |       -5      +0.47         +3 | 1         +3 |       1 | uniform(-5,3) 
#   7|D  | k5         |       -5       -2.5         +3 | 1    +0.0031 |       1 | uniform(-5,3) 
#   8|D  | kb1        |       -5       -1.8         +3 | 1     +0.017 |       1 | uniform(-5,3) 
#   9|D  | kb2        |       -5       -2.2         +3 | 1    +0.0058 |       1 | uniform(-5,3) 
#  10|D  | kc1        |       -5       -2.8         +3 | 1    +0.0017 |       1 | uniform(-5,3) 
     |   |            |                                |              |         |      
#  11|D  | kc2        |       -5       -2.2         +3 | 1     +0.007 |       1 | uniform(-5,3) 
#  12|D  | kc4        |       -5       -2.2         +3 | 1    +0.0061 |       1 | uniform(-5,3) 
#  13|D  | szea       |       -5      -0.29         +3 | 1      +0.52 |       1 | uniform(-5,3) 


Please write down a short documentation:

The model describes cleavage of beta-carotines in Arabidopsis.

The publisehd model was set up in D2D by Clemens Kreutz together with Mirjam Fehling-Kaschek (MFK), one of the coauthors of the above publication.
Experimental errors were calculated by MFK by estimating an error model to standard deviations determined from experimetnal replicates.


