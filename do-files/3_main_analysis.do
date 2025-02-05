/*==============================================================================
DO FILE NAME:			3_analysis.do
DATE: 					04/02/2025
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	runs main analysis

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\3_analysis_$S_DATE.log", append



use "$savedir\an_dm_main_analysis", clear
* Model treatment allocation on the set of confounding variables  
logistic antivegf age_at_index gender eth5 imd smokstatus bmi_cat yrs_dm dm_type drug_dm_count yrs_retinopathy yrs_eye_dis bl_amputation bl_neuropathy bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_arni index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior

* Estimate propensity scores
predict propensity

graph tw kdensity propensity if t == 0 || kdensity propensity if t == 1
estat gof, group(10) table

* Calculate ATT weights 
gen att_weight = antivegf + (1-antivegf)*(propensity/(1-propensity))

* Fit weighted Cox regression w/ robust standard errors
stset end_fu [pweight=att_weight], failure(egfr_40) origin(index_date) enter(index_date) scale(365.25) id(patid)
stcox antivegf, vce(robust)
