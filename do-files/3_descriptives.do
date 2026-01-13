/*==============================================================================
DO FILE NAME:			2_descriptives.do
DATE: 					20/11/2024
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	generates descriptive statistics

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\2_descriptives_$S_DATE.log", append


* Cohort characteristics
* Diabetes groups
use "$savedir\an_dm_main_analysis", clear
drop if out==1
dtable fu_time age_at_index i.gender i.eth5 i.bmi_cat i.smokstatus i.imd yrs_dm i.dm_type_f i.index_mtf i.drug_t2dm i.index_ins i.drug_aaa i.index_mra i.index_antiplatelet i.drug_antihyp i.index_oac i.index_ppi i.index_statin i.index_nsaid i.bl_retinopathy i.bl_eye_dis yrs_retinopathy i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_ncm i.bl_neph_syndrome i.bl_rrt i.bl_stroke i.bl_zoster i.bl_hf i.bl_pad i.bl_copd i.bl_neuropathy i.bl_amputation  prior_egfr i.prior_egfr_cat i.prior_ckd i.prior_acr_cat tot_appts_yr_prior scr_meas_yr_prior, by(antivegf) export($projdir\output\dm_table_1.xlsx, replace)

* Number of outcomes
dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_zoster, by(antivegf) export($projdir\output\dm_outcomes_desc.xlsx, replace)


* No diabetes groups
use "$savedir\an_nodm_main_analysis", clear
drop if out==1
gen fu_time_yrs  = fu_time/365.25
dtable fu_time age_at_index i.gender i.eth5 i.bmi_cat i.smokstatus i.imd i.drug_aaa i.index_mra i.index_antiplatelet i.drug_antihyp i.index_oac i.index_ppi i.index_statin i.index_nsaid i.bl_retinopathy i.bl_eye_dis i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_ncm i.bl_neph_syndrome i.bl_rrt i.bl_stroke i.bl_zoster i.bl_hf i.bl_pad i.bl_copd i.bl_neuropathy i.bl_amputation  prior_egfr i.prior_egfr_cat i.prior_ckd i.prior_acr_cat tot_appts_yr_prior scr_meas_yr_prior, by(antivegf) export($projdir\output\no_dm_table_1.xlsx, replace)

* Number of outcomes
dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_zoster, by(antivegf) export($projdir\output\nodm_outcomes_desc.xlsx, replace)

* Number of antivegf injections 
* Need to use only people in kept in the analysis
use "$savedir\an_dm_main_analysis", clear
append using "$savedir\an_nodm_main_analysis"
keep if antivegf==1
keep patid out bl_dm
drop if out==1
merge 1:m patid using "$savedir\antivegf\injections", keep(match) keepusing(yr_since_fu number_antivegf_yr fu_yr)
dtable number_antivegf_yr if bl_dm==0, by(yr_since_fu) continuous(number_antivegf_yr, statistics(median q1 q3 min max)) export($projdir\output\no_dm_table_2.xlsx, replace)
dtable number_antivegf_yr if bl_dm==1, by(yr_since_fu) continuous(number_antivegf_yr, statistics(median q1 q3 min max)) export($projdir\output\dm_table_2.xlsx, replace)
bys patid: gen n=_n
tab fu_yr bl_dm if n==1


* Number of outcomes within strata
foreach group in dm nodm {
	use "$savedir\an_`group'_main_analysis", clear
	drop if out==1
* eGFR
	foreach strata in prior_egfr_cat prior_acr_cat eth5 imd gender {
		dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_zoster if antivegf==1, by(`strata') export($projdir\output\\`group'_outcomes_desc_exposed_`strata'.xlsx, replace)
		dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_zoster if antivegf==0, by(`strata') export($projdir\output\\`group'_outcomes_desc_unexposed_`strata'.xlsx, replace)
	}
}




