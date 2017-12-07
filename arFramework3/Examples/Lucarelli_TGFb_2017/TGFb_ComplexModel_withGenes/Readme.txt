Original folder name: /dsk/home/ck96/2013_TGFb/19_geneExpression_NachMS 
Date: 13-Jun-2017 


Model states and right hand side of the ODEs:
TGFb		-Rec*Rec_act*TGFb
Rec		-Rec*Rec_act*TGFb
TGFb_pRec		Rec*Rec_act*TGFb - TGFb_pRec*pRec_degind
S2		- S2*k_on_u*S4^2 + S2_S4_S4*kdiss_SS + S_dephos*pS2 - S2*S_phos*TGFb_pRec
S3		S_dephos*pS3 - S3*S_phos*TGFb_pRec
S4		- 3*khomo4*S4^3 - 2*k_244*S4^2*ppS2 - 2*k_344*S4^2*ppS3 - 2*S2*k_on_u*S4^2 - k_224*S4*ppS2^2 - k_234*S4*ppS2*ppS3 - k_334*S4*ppS3^2 + 2*S2_S4_S4*kdiss_SS + 3*S4_S4_S4*kdiss_SS + 2*S_dephosphos*ppS2_S4_S4 + 2*S_dephosphos*ppS3_S4_S4 + 2*S_dephosphos*ppS2_ppS2_S4 + 2*S_dephosphos*ppS2_ppS3_S4 + 2*S_dephosphos*ppS3_ppS3_S4
S2_S4_S4		S2*k_on_u*S4^2 - S2_S4_S4*kdiss_SS
ppS2_ppS2_ppS2		khomo2*ppS2^3 - 3*S_dephosphos*ppS2_ppS2_ppS2
ppS3_ppS3_ppS3		khomo3*ppS3^3 - 3*S_dephosphos*ppS3_ppS3_ppS3
S4_S4_S4		khomo4*S4^3 - S4_S4_S4*kdiss_SS
pS2		S_dephosphos*ppS2 - S_dephos*pS2 + S_dephosphos*ppS2_S4_S4 + 2*S_dephosphos*ppS2_ppS2_S4 + S_dephosphos*ppS2_ppS3_S4 + 3*S_dephosphos*ppS2_ppS2_ppS2 + 2*S_dephosphos*ppS2_ppS2_ppS3 + S_dephosphos*ppS2_ppS3_ppS3
pS3		S_dephosphos*ppS3 - S_dephos*pS3 + S_dephosphos*ppS3_S4_S4 + S_dephosphos*ppS2_ppS3_S4 + 2*S_dephosphos*ppS3_ppS3_S4 + S_dephosphos*ppS2_ppS2_ppS3 + 2*S_dephosphos*ppS2_ppS3_ppS3 + 3*S_dephosphos*ppS3_ppS3_ppS3
ppS2		- k_244*S4^2*ppS2 - 2*k_224*S4*ppS2^2 - k_234*S4*ppS2*ppS3 - 3*khomo2*ppS2^3 - 2*k_223*ppS2^2*ppS3 - k_233*ppS2*ppS3^2 - S_dephosphos*ppS2 + 2*S_dephosphos*ppS2_ppS2_S4 + S_dephosphos*ppS2_ppS3_S4 + 6*S_dephosphos*ppS2_ppS2_ppS2 + 4*S_dephosphos*ppS2_ppS2_ppS3 + 2*S_dephosphos*ppS2_ppS3_ppS3 + S2*S_phos*TGFb_pRec
ppS3		- k_344*S4^2*ppS3 - k_234*S4*ppS2*ppS3 - 2*k_334*S4*ppS3^2 - k_223*ppS2^2*ppS3 - 2*k_233*ppS2*ppS3^2 - 3*khomo3*ppS3^3 - S_dephosphos*ppS3 + S_dephosphos*ppS2_ppS3_S4 + 2*S_dephosphos*ppS3_ppS3_S4 + 2*S_dephosphos*ppS2_ppS2_ppS3 + 4*S_dephosphos*ppS2_ppS3_ppS3 + 6*S_dephosphos*ppS3_ppS3_ppS3 + S3*S_phos*TGFb_pRec
ppS2_ppS2_S4		S4*k_224*ppS2^2 - 2*S_dephosphos*ppS2_ppS2_S4
ppS2_ppS2_ppS3		k_223*ppS3*ppS2^2 - 3*S_dephosphos*ppS2_ppS2_ppS3
ppS2_ppS3_ppS3		k_233*ppS2*ppS3^2 - 3*S_dephosphos*ppS2_ppS3_ppS3
ppS3_ppS3_S4		S4*k_334*ppS3^2 - 2*S_dephosphos*ppS3_ppS3_S4
ppS2_ppS3_S4		S4*k_234*ppS2*ppS3 - 2*S_dephosphos*ppS2_ppS3_S4
ppS3_S4_S4		k_344*ppS3*S4^2 - S_dephosphos*ppS3_S4_S4
ppS2_S4_S4		k_244*ppS2*S4^2 - S_dephosphos*ppS2_S4_S4
geneA		(geneA_turn + geneA_act2*ppS2_S4_S4 + geneA_act1*ppS2_ppS3_S4 + geneA_act3*ppS2_ppS3_ppS3)/(geneA_inh2*ppS2_S4_S4 + geneA_inh1*ppS2_ppS3_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1) - geneA*geneA_turn
geneB		(geneB_turn + geneB_act2*ppS2_S4_S4 + geneB_act1*ppS2_ppS3_S4 + geneB_act3*ppS2_ppS3_ppS3)/(geneB_inh2*ppS2_S4_S4 + geneB_inh1*ppS2_ppS3_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1) - geneB*geneB_turn
geneC		(geneC_turn + geneC_act2*ppS2_S4_S4 + geneC_act1*ppS2_ppS3_S4 + geneC_act3*ppS2_ppS3_ppS3)/(geneC_inh2*ppS2_S4_S4 + geneC_inh1*ppS2_ppS3_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1) - geneC*geneC_turn
geneD		(geneD_turn + geneD_act2*ppS2_S4_S4 + geneD_act1*ppS2_ppS3_S4 + geneD_act3*ppS2_ppS3_ppS3)/(geneD_inh2*ppS2_S4_S4 + geneD_inh1*ppS2_ppS3_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1) - geneD*geneD_turn
geneE		(geneE_turn + geneE_act2*ppS2_S4_S4 + geneE_act1*ppS2_ppS3_S4 + geneE_act3*ppS2_ppS3_ppS3)/(geneE_inh2*ppS2_S4_S4 + geneE_inh1*ppS2_ppS3_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1) - geneE*geneE_turn
geneF		(geneF_turn + geneF_act2*ppS2_S4_S4 + geneF_act1*ppS2_ppS3_S4 + geneF_act3*ppS2_ppS3_ppS3)/(geneF_inh2*ppS2_S4_S4 + geneF_inh1*ppS2_ppS3_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1) - geneF*geneF_turn
geneG		(geneG_turn + geneG_act2*ppS2_S4_S4 + geneG_act1*ppS2_ppS3_S4 + geneG_act3*ppS2_ppS3_ppS3)/(geneG_inh2*ppS2_S4_S4 + geneG_inh1*ppS2_ppS3_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1) - geneG*geneG_turn
geneH		(geneH_turn + geneH_act2*ppS2_S4_S4 + geneH_act1*ppS2_ppS3_S4 + geneH_act3*ppS2_ppS3_ppS3)/(geneH_inh2*ppS2_S4_S4 + geneH_inh1*ppS2_ppS3_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1) - geneH*geneH_turn
geneI		(geneI_turn + geneI_act2*ppS2_S4_S4 + geneI_act1*ppS2_ppS3_S4 + geneI_act3*ppS2_ppS3_ppS3)/(geneI_inh2*ppS2_S4_S4 + geneI_inh1*ppS2_ppS3_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1) - geneI*geneI_turn
geneJ		(geneJ_turn + geneJ_act2*ppS2_S4_S4 + geneJ_act1*ppS2_ppS3_S4 + geneJ_act3*ppS2_ppS3_ppS3)/(geneJ_inh2*ppS2_S4_S4 + geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1) - geneJ*geneJ_turn
geneK		(geneK_turn + geneK_act2*ppS2_S4_S4 + geneK_act1*ppS2_ppS3_S4 + geneK_act3*ppS2_ppS3_ppS3)/(geneK_inh2*ppS2_S4_S4 + geneK_inh1*ppS2_ppS3_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1) - geneK*geneK_turn
geneL		(geneL_turn + geneL_act2*ppS2_S4_S4 + geneL_act1*ppS2_ppS3_S4 + geneL_act3*ppS2_ppS3_ppS3)/(geneL_inh2*ppS2_S4_S4 + geneL_inh1*ppS2_ppS3_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1) - geneL*geneL_turn


Observables and Errors:
S2IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS4:		fy= ( ( ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS4:		fy= ( ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS2:		fy= ( ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS3:		fy= ( ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS4:		fy= ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 +  ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS4:		fy= ( ( ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_pS2:		fy= ( ( pS2 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_S2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS2:		fy= ( ( pS2 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_S3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_pS3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_ppS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3) + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS4:		fy= ( ( ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS4:		fy= ( ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS2:		fy= ( ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS3:		fy= ( ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS4:		fy= ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 +  ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS4:		fy= ( ( ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_pS2:		fy= ( ( pS2 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_S2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS2:		fy= ( ( pS2 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_S3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_pS3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_ppS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3) + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS4:		fy= ( ( ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS4:		fy= ( ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS2:		fy= ( ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS3:		fy= ( ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS4:		fy= ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS2:		fy= ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 +  ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS3:		fy= ( ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S23IP_tS4:		fy= ( ( ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) / ( ( S2 +  S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_pS2:		fy= ( ( pS2 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_S2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS2:		fy= ( ( 0 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S2:		fy= ( ( S2 + S2_S4_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS2:		fy= ( ( pS2 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS2:		fy= ( (( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( S2 + S2_S4_S4 + pS2 +( ppS2 + 3 * ppS2_ppS2_ppS2 )+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_S3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_pS3:		fy= ( ( 0 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S2IP_ppS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S3IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S3:		fy= ( ( S3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS3:		fy= ( ( pS3 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3) + ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S2:		fy= ( ( S2  +  S2_S4_S4 ) / ( S2  +  S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_pS2:		fy= ( ( pS2 ) / ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS2:		fy= ( ((ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( S2  +  S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_S3:		fy= ( ( S3 ) / ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) * 100		fystd=5
S23IP_pS3:		fy= ( ( pS3 ) / ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_ppS3:		fy= ( ((ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) )  * 100		fystd=5
S23IP_tS2:		fy=(S2 + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4) 		fystd=15
S4IP_tS4:		fy=((S4 + 3*S4_S4_S4) + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4  + ppS2_ppS3_S4)		fystd=15
S23IP_tS2:		fy=(S2 + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2) + 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4)/((S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4)+(S2 + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4)) * 100 		fystd=1
S23IP_tS3:		fy=(S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4)/((S3 + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4)+(S2 + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4)) * 100 		fystd=1
TGFbR:		fy=(Rec + TGFb_pRec)		fystd=0.092
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
Ski:		fy=1e-6 + geneA		fystd=sd_Ski_nExpID100
Skil:		fy=1e-6 + geneB		fystd=sd_Skil_nExpID100
Dnmt3a:		fy=1e-6 + geneC		fystd=sd_Dnmt3a_nExpID100
Sox4:		fy=1e-6 + geneD		fystd=sd_Sox4_nExpID100
Jun:		fy=1e-6 + geneE		fystd=sd_Jun_nExpID100
Smad7:		fy=1e-6 + geneF		fystd=sd_Smad7_nExpID100
Klf10:		fy=1e-6 + geneG		fystd=sd_Klf10_nExpID100
Bmp4:		fy=1e-6 + geneH		fystd=sd_Bmp4_nExpID100
Cxcl15:		fy=1e-6 + geneI		fystd=sd_Cxcl15_nExpID100
Dusp5:		fy=1e-6 + geneJ		fystd=sd_Dusp5_nExpID100
Tgfa:		fy=1e-6 + geneK		fystd=sd_Tgfa_nExpID100
Pdk4:		fy=1e-6 + geneL		fystd=sd_Pdk4_nExpID100
S2IP_tS2:		fy= ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS3:		fy= ( ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S2  +  S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_tS4:		fy= ( ( ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( 2 * S2_S4_S4 + ppS2_ppS2_S4 + 2 * ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS2:		fy= ( ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 ) / ( ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS3:		fy= ( ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S3IP_tS4:		fy= ( ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( S3  + pS3 +(ppS3+3*ppS3_ppS3_ppS3)+ ppS2_ppS2_ppS3 + 2 * ppS2_ppS3_ppS3 + 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS2_ppS2_ppS3 + ppS2_ppS3_ppS3 + ppS2_ppS3_S4 )  +  ( ppS3_ppS3_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS2:		fy= ( ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS3:		fy= ( ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S4IP_tS4:		fy= ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 ) / ( ( ( S4 + 3 * S4_S4_S4 )  + 2 * S2_S4_S4 + ppS2_ppS2_S4 + ppS3_ppS3_S4 + 2 * ppS2_S4_S4 + 2 * ppS3_S4_S4 + ppS2_ppS3_S4 )  +  ( S2_S4_S4 + 2 * ppS2_ppS2_S4 + ppS2_S4_S4 + ppS2_ppS3_S4 )  +  ( 2 * ppS3_ppS3_S4 + ppS3_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_S2:		fy= ( ( S2  + S2_S4_S4 ) / ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_pS2:		fy= ( ( pS2 ) / ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5
S2IP_ppS2:		fy= ( ((ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) / ( ( S2  + S2_S4_S4 + pS2 +(ppS2+3*ppS2_ppS2_ppS2)+ 2 * ppS2_ppS2_ppS3 + 2 * ppS2_ppS2_S4 + ppS2_ppS3_ppS3 + ppS2_S4_S4 + ppS2_ppS3_S4 ) ) )  * 100		fystd=5


Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model

           name                  lb       value       ub          10^value        fitted   prior
#   1|D  | Rec_act             |       -5       -2.2         +3 | 1     +0.006 |       1 | uniform(-5,3) 
#   2|DI | S2tot               |       -5       +2.2         +3 | 1   +1.4e+02 |       1 | uniform(-5,3) 
#   3|DI | S3tot               |       -5       +1.2         +3 | 1        +16 |       1 | uniform(-5,3) 
#   4|DI | S4tot               |       -5       +1.8         +3 | 1        +67 |       1 | uniform(-5,3) 
#   5|D  | S_dephos            |       -5      -0.54         +3 | 1      +0.29 |       1 | uniform(-5,3) 
#   6|D  | S_dephosphos        |       -5       -1.3         +3 | 1     +0.049 |       1 | uniform(-5,3) 
#   7|D  | S_phos              |       -5      -0.42         +3 | 1      +0.38 |       1 | uniform(-5,3) 
#   8|D  | geneA_act1          |       -5       -1.9         +3 | 1     +0.014 |       1 | uniform(-5,3) 
#   9|D  | geneA_act2          |       -5       -3.1         +3 | 1   +0.00081 |       1 | uniform(-5,3) 
#  10|D  | geneA_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  11|D  | geneA_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  12|D  | geneA_inh2          |       -5       -1.3         +3 | 1     +0.047 |       1 | uniform(-5,3) 
#  13|D  | geneA_inh3          |       -5       -1.6         +3 | 1     +0.027 |       1 | uniform(-5,3) 
#  14|D  | geneA_turn          |       -5       -2.3         +3 | 1    +0.0046 |       1 | uniform(-5,3) 
#  15|D  | geneB_act1          |       -5      -0.33         +3 | 1      +0.46 |       1 | uniform(-5,3) 
#  16|D  | geneB_act2          |       -5       -1.1         +3 | 1     +0.075 |       1 | uniform(-5,3) 
#  17|D  | geneB_act3          |       -5       -1.5         +3 | 1     +0.031 |       1 | uniform(-5,3) 
#  18|D  | geneB_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  19|D  | geneB_inh2          |       -5      -0.12         +3 | 1      +0.75 |       1 | uniform(-5,3) 
#  20|D  | geneB_inh3          |       -5      -0.39         +3 | 1      +0.41 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  21|D  | geneB_turn          |       -5       -1.9         +3 | 1     +0.011 |       1 | uniform(-5,3) 
#  22|D  | geneC_act1          |       -5       -2.1         +3 | 1    +0.0072 |       1 | uniform(-5,3) 
#  23|D  | geneC_act2          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  24|D  | geneC_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  25|D  | geneC_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  26|D  | geneC_inh2          |       -5       -1.7         +3 | 1     +0.019 |       1 | uniform(-5,3) 
#  27|D  | geneC_inh3          |       -5       -1.4         +3 | 1     +0.036 |       1 | uniform(-5,3) 
#  28|D  | geneC_turn          |       -5       -2.4         +3 | 1    +0.0038 |       1 | uniform(-5,3) 
#  29|D  | geneD_act1          |       -5         -2         +3 | 1     +0.011 |       1 | uniform(-5,3) 
#  30|D  | geneD_act2          |       -5       -3.6         +3 | 1   +0.00027 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  31|D  | geneD_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  32|D  | geneD_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  33|D  | geneD_inh2          |       -5       -1.1         +3 | 1     +0.087 |       1 | uniform(-5,3) 
#  34|D  | geneD_inh3          |       -5      -0.99         +3 | 1       +0.1 |       1 | uniform(-5,3) 
#  35|D  | geneD_turn          |       -5       -3.1         +3 | 1    +0.0008 |       1 | uniform(-5,3) 
#  36|D  | geneE_act1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  37|D  | geneE_act2          |       -5    +0.0083         +3 | 1         +1 |       1 | uniform(-5,3) 
#  38|D  | geneE_act3          |       -5      +0.78         +3 | 1         +6 |       1 | uniform(-5,3) 
#  39|D  | geneE_inh1          |       -5      +0.91         +3 | 1       +8.2 |       1 | uniform(-5,3) 
#  40|D  | geneE_inh2          |       -5      +0.16         +3 | 1       +1.5 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  41|D  | geneE_inh3          |       -5      +0.99         +3 | 1       +9.8 |       1 | uniform(-5,3) 
#  42|D  | geneE_turn          |       -5      -0.92         +3 | 1      +0.12 |       1 | uniform(-5,3) 
#  43|D  | geneF_act1          |       -5         +3         +3 | 1     +1e+03 |       1 | uniform(-5,3) 
#  44|D  | geneF_act2          |       -5      -0.67         +3 | 1      +0.21 |       1 | uniform(-5,3) 
#  45|D  | geneF_act3          |       -5       +1.3         +3 | 1        +21 |       1 | uniform(-5,3) 
#  46|D  | geneF_inh1          |       -5      +0.56         +3 | 1       +3.7 |       1 | uniform(-5,3) 
#  47|D  | geneF_inh2          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  48|D  | geneF_inh3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  49|D  | geneF_turn          |       -5       +1.6         +3 | 1        +36 |       1 | uniform(-5,3) 
#  50|D  | geneG_act1          |       -5         +3         +3 | 1     +1e+03 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  51|D  | geneG_act2          |       -5         +2         +3 | 1        +91 |       1 | uniform(-5,3) 
#  52|D  | geneG_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  53|D  | geneG_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  54|D  | geneG_inh2          |       -5       -1.7         +3 | 1     +0.018 |       1 | uniform(-5,3) 
#  55|D  | geneG_inh3          |       -5       -3.1         +3 | 1   +0.00076 |       1 | uniform(-5,3) 
#  56|D  | geneG_turn          |       -5       +2.7         +3 | 1   +4.7e+02 |       1 | uniform(-5,3) 
#  57|D  | geneH_act1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  58|D  | geneH_act2          |       -5         +2         +3 | 1        +90 |       1 | uniform(-5,3) 
#  59|D  | geneH_act3          |       -5         +3         +3 | 1     +1e+03 |       1 | uniform(-5,3) 
#  60|D  | geneH_inh1          |       -5         +3         +3 | 1     +1e+03 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  61|D  | geneH_inh2          |       -5       +1.3         +3 | 1        +19 |       1 | uniform(-5,3) 
#  62|D  | geneH_inh3          |       -5       +2.3         +3 | 1   +2.2e+02 |       1 | uniform(-5,3) 
#  63|D  | geneH_turn          |       -5    -0.0037         +3 | 1      +0.99 |       1 | uniform(-5,3) 
#  64|D  | geneI_act1          |       -5      +0.97         +3 | 1       +9.2 |       1 | uniform(-5,3) 
#  65|D  | geneI_act2          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  66|D  | geneI_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  67|D  | geneI_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  68|D  | geneI_inh2          |       -5       +2.5         +3 | 1   +3.2e+02 |       1 | uniform(-5,3) 
#  69|D  | geneI_inh3          |       -5         +3         +3 | 1     +1e+03 |       1 | uniform(-5,3) 
#  70|D  | geneI_turn          |       -5       -2.5         +3 | 1    +0.0033 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  71|D  | geneJ_act1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  72|D  | geneJ_act2          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  73|D  | geneJ_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  74|D  | geneJ_inh1          |       -5      -0.33         +3 | 1      +0.47 |       1 | uniform(-5,3) 
#  75|D  | geneJ_inh2          |       -5       -1.9         +3 | 1     +0.012 |       1 | uniform(-5,3) 
#  76|D  | geneJ_inh3          |       -5         -2         +3 | 1     +0.011 |       1 | uniform(-5,3) 
#  77|D  | geneJ_turn          |       -5       -1.8         +3 | 1     +0.014 |       1 | uniform(-5,3) 
#  78|D  | geneK_act1          |       -5       -1.9         +3 | 1     +0.013 |       1 | uniform(-5,3) 
#  79|D  | geneK_act2          |       -5       -3.5         +3 | 1   +0.00029 |       1 | uniform(-5,3) 
#  80|D  | geneK_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  81|D  | geneK_inh1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  82|D  | geneK_inh2          |       -5      +0.14         +3 | 1       +1.4 |       1 | uniform(-5,3) 
#  83|D  | geneK_inh3          |       -5       +0.4         +3 | 1       +2.5 |       1 | uniform(-5,3) 
#  84|D  | geneK_turn          |       -5       -2.4         +3 | 1    +0.0036 |       1 | uniform(-5,3) 
#  85|D  | geneL_act1          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  86|D  | geneL_act2          |       -5       -3.3         +3 | 1   +0.00049 |       1 | uniform(-5,3) 
#  87|D  | geneL_act3          |       -5         -5         +3 | 1     +1e-05 |       1 | uniform(-5,3) 
#  88|D  | geneL_inh1          |       -5      +0.11         +3 | 1       +1.3 |       1 | uniform(-5,3) 
#  89|D  | geneL_inh2          |       -5       -1.2         +3 | 1     +0.064 |       1 | uniform(-5,3) 
#  90|D  | geneL_inh3          |       -5      -0.85         +3 | 1      +0.14 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
#  91|D  | geneL_turn          |       -5       -1.8         +3 | 1     +0.015 |       1 | uniform(-5,3) 
#  92|DI | init_Rec            |       -5      +0.26         +3 | 1       +1.8 |       1 | uniform(-5,3) 
#  93|D  | k_223               |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
#  94|D  | k_224               |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
#  95|D  | k_233               |       -8      -0.83         +3 | 1      +0.15 |       1 | uniform(-8,3) 
#  96|D  | k_234               |       -8       -3.4         +3 | 1    +0.0004 |       1 | uniform(-8,3) 
#  97|D  | k_244               |       -8       -5.1         +3 | 1   +8.6e-06 |       1 | uniform(-8,3) 
#  98|D  | k_334               |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
#  99|D  | k_344               |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
# 100|D  | k_on_u              |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
     |   |                     |                                |              |         |      
# 101|D  | kdiss_SS            |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
# 102|D  | khomo2              |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
# 103|D  | khomo3              |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
# 104|D  | khomo4              |       +0         +0         +3 | 0         +0 |       0 | uniform(0,3) 
# 105|D  | pRec_degind         |       -5       -1.4         +3 | 1      +0.04 |       1 | uniform(-5,3) 
# 106|  E| sd_Bmp4_nExpID100   |       -5      -0.85         +3 | 1      +0.14 |       1 | uniform(-5,3) 
# 107|  E| sd_Cxcl15_nExpID100 |       -5      -0.87         +3 | 1      +0.13 |       1 | uniform(-5,3) 
# 108|  E| sd_Dnmt3a_nExpID100 |       -5         -1         +3 | 1     +0.097 |       1 | uniform(-5,3) 
# 109|  E| sd_Dusp5_nExpID100  |       -5      -0.98         +3 | 1      +0.11 |       1 | uniform(-5,3) 
# 110|  E| sd_Jun_nExpID100    |       -5      -0.78         +3 | 1      +0.16 |       1 | uniform(-5,3) 
     |   |                     |                                |              |         |      
# 111|  E| sd_Klf10_nExpID100  |       -5      -0.93         +3 | 1      +0.12 |       1 | uniform(-5,3) 
# 112|  E| sd_Pdk4_nExpID100   |       -5       -1.1         +3 | 1     +0.074 |       1 | uniform(-5,3) 
# 113|  E| sd_Ski_nExpID100    |       -5      -0.99         +3 | 1       +0.1 |       1 | uniform(-5,3) 
# 114|  E| sd_Skil_nExpID100   |       -5      -0.76         +3 | 1      +0.17 |       1 | uniform(-5,3) 
# 115|  E| sd_Smad7_nExpID100  |       -5       -1.2         +3 | 1      +0.07 |       1 | uniform(-5,3) 
# 116|  E| sd_Sox4_nExpID100   |       -5       -0.9         +3 | 1      +0.13 |       1 | uniform(-5,3) 
# 117|  E| sd_Tgfa_nExpID100   |       -5       -1.1         +3 | 1     +0.077 |       1 | uniform(-5,3) 


Please write down a short documentation:

TGFb complex model after selecting the relevant complexes and after linking
them to genes.

Original workspace: /dsk/home/ck96/2013_TGFb/19_geneExpression_NachMS/Results/20150324T093446_5_Recompile      
( -2*log(L)=213.833  N=1755  #p=117  #fitted=108 ./Results/ #prior=  0  errors fitted              # old PLE=108 )

