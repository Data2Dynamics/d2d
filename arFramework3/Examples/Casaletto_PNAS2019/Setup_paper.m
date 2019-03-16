
arInit;

ar.config.checkForNegFluxes = false;

arLoadModel('dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('egfr_cmet_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('mm151_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('mm151_egfr_cMet_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('no_ta_cMet_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy2_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy3_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy3_ext_nm_dr','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy4_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy5_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy6_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy7_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy8_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('dummy9_ext_nm','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('jnj_bs_nm_dr_snu5','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_egfr_nm_dr_snu5','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_nm_dr_snu5','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_egfr_nm_dr_snu5','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('jnj_bs_nm_dr_h1993','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_egfr_nm_dr_h1993','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_nm_dr_h1993','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_egfr_nm_dr_h1993','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('jnj_bs_nm_dr_h292','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_egfr_nm_dr_h292','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_nm_dr_h292','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jnj_met_egfr_nm_dr_h292','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('ssu_data','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('ssu_data_mm131','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('jc_data_mm131_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('jc_data_metmab_nm','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('jc_data_metmab_nm_h441','dimer_model_paper_make_clean_koff_no_ratios_der');
% arLoadData('jc_data_metmab_nm_h747','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('jc_data_mm131_nm_ID','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jc_data_metmab_nm_ID','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jc_data_metmab_nm_h441_ID','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('jc_data_metmab_nm_h747_ID','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('hgf_egf_combo_inhibition_2_AR1','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('hgf_egf_combo_inhibition_2_AR2','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('hgf_egf_combo_inhibition_2_AR3','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('hgf_egf_combo_inhibition_2_AR4','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('hgf_egf_combo_inhibition_2_AR5','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('egfr_cmet_nm_dummy','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('figure_8_simu_1','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_8_simu_2','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_8_simu_3','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('figure_9_simu_1','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_2','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_3','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_4','dimer_model_paper_make_clean_koff_no_ratios_der');

arLoadData('figure_9_simu_1_DR','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_2_DR','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_3_DR','dimer_model_paper_make_clean_koff_no_ratios_der');
arLoadData('figure_9_simu_4_DR','dimer_model_paper_make_clean_koff_no_ratios_der');

arCompileAll;

arMergePlotMulti(1,{'mm151_egfr_cMet_nm','egfr_cmet_nm','no_ta_cMet_nm','mm151_nm'},{'mm151_egfr_cMet_nm','egfr_cMet_nm','no_ta_cMet_nm','mm151_nm'},'combi_condi')

arMergePlotMulti(1,{'jc_data_metmab_nm_h747_ID_ID1','jc_data_metmab_nm_h747_ID_ID2','jc_data_metmab_nm_h747_ID_ID3','jc_data_metmab_nm_h747_ID_ID4','jc_data_metmab_nm_h747_ID_ID5','jc_data_metmab_nm_h747_ID_ID6'},{'jc_data_metmab_nm_h747_ID_ID1','jc_data_metmab_nm_h747_ID_ID2','jc_data_metmab_nm_h747_ID_ID3','jc_data_metmab_nm_h747_ID_ID4','jc_data_metmab_nm_h747_ID_ID5','jc_data_metmab_nm_h747_ID_ID6'},'metmab_nm_h747')
arMergePlotMulti(1,{'jc_data_metmab_nm_h441_ID_ID1','jc_data_metmab_nm_h441_ID_ID2','jc_data_metmab_nm_h441_ID_ID3','jc_data_metmab_nm_h441_ID_ID4'},{'jc_data_metmab_nm_h441_ID_ID1','jc_data_metmab_nm_h441_ID_ID2','jc_data_metmab_nm_h441_ID_ID3','jc_data_metmab_nm_h441_ID_ID4'},'metmab_nm_h441')
arMergePlotMulti(1,{'jc_data_metmab_nm_ID_ID1','jc_data_metmab_nm_ID_ID2','jc_data_metmab_nm_ID_ID3','jc_data_metmab_nm_ID_ID4'},{'jc_data_metmab_nm_ID_ID1','jc_data_metmab_nm_ID_ID2','jc_data_metmab_nm_ID_ID3','jc_data_metmab_nm_ID_ID4'},'metmab_nm_a549')
arMergePlotMulti(1,{'jc_data_mm131_nm_ID_ID1','jc_data_mm131_nm_ID_ID2','jc_data_mm131_nm_ID_ID3','jc_data_mm131_nm_ID_ID4'},{'jc_data_mm131_nm_ID_ID1','jc_data_mm131_nm_ID_ID2','jc_data_mm131_nm_ID_ID3','jc_data_mm131_nm_ID_ID4'},'mm131_jc')

arMergePlotMulti(1,{'hgf_egf_combo_inhibition_2_AR1','hgf_egf_combo_inhibition_2_AR2','hgf_egf_combo_inhibition_2_AR3','hgf_egf_combo_inhibition_2_AR4','hgf_egf_combo_inhibition_2_AR5'},{'MM151','MetmAB','MM131','MM151 + MetmAB','MM151 + MM131'},'ssu_hgf_egf_combi')

% arMergePlotMulti(1,{'jnj_bs_nm_dr_snu5','jnj_egfr_nm_dr_snu5','jnj_met_nm_dr_snu5','jnj_met_egfr_nm_dr_snu5'},{'BsAb','EGFR_blind','Met_blind','EGFR_blind + Met_blind'},'JNJ_combi_snu5')
% arMergePlotMulti(1,{'jnj_bs_nm_dr_h1993','jnj_egfr_nm_dr_h1993','jnj_met_nm_dr_h1993','jnj_met_egfr_nm_dr_h1993'},{'BsAb','EGFR_blind','Met_blind','EGFR_blind + Met_blind'},'JNJ_combi_h1993')
% arMergePlotMulti(1,{'jnj_bs_nm_dr_h292','jnj_egfr_nm_dr_h292','jnj_met_nm_dr_h292','jnj_met_egfr_nm_dr_h292'},{'BsAb','EGFR_blind','Met_blind','EGFR_blind + Met_blind'},'JNJ_combi_h292')

arMergePlotMulti(1,{'jnj_bs_nm_dr_snu5','jnj_met_nm_dr_snu5'},{'BsAb','Met_blind'},'JNJ_combi_snu5')
arMergePlotMulti(1,{'jnj_bs_nm_dr_h1993','jnj_met_nm_dr_h1993'},{'BsAb','Met_blind'},'JNJ_combi_h1993')
arMergePlotMulti(1,{'jnj_bs_nm_dr_h292','jnj_met_nm_dr_h292'},{'BsAb','Met_blind'},'JNJ_combi_h292')


arMergePlotMulti(1,{'figure_8_simu_1','figure_8_simu_2','figure_8_simu_3'},{'EGFR_Met_BS','MM151 + Metmab','MM151 + MM131'},'figure_8_simu')
arMergePlotMulti(1,{'figure_9_simu_1','figure_9_simu_2','figure_9_simu_3','figure_9_simu_4'},{'Control','MM151','MM131','MM151 + MM131'},'figure_9_simu')
arMergePlotMulti(1,{'figure_9_simu_1_DR','figure_9_simu_2_DR','figure_9_simu_3_DR','figure_9_simu_4_DR'},{'Control','MM151','MM131','MM151 + MM131'},'figure_9_simu_DR')
