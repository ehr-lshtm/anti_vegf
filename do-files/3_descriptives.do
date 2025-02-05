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

use "$savedir\antivegf\cr_vars", clear
keep if bl_dm==1
gen antivegf=1
append using "$savedir\photocoag\cr_vars"
replace antivegf = 0 if antivegf==.

dtable age_at_index i.gender i.eth5 i.bmi_cat i.smokstatus yrs_dm i.dm_type_f i.index_acarbose-index_nsaid i.bl_retinopathy i.bl_eye_dis i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_ncm i.bl_neph_syndrome i.bl_rrt i.bl_stroke i.bl_zoster i.bl_hf i.bl_pad i.bl_copd i.bl_neuropathy i.bl_amputation  prior_egfr i.prior_egfr_cat i.prior_ckd i.prior_acr_cat tot_appts_yr_prior scr_meas_yr_prior, by(antivegf) export($projdir\output\dm_table_1.xlsx, replace)

* No diabetes groups

use "$savedir\antivegf\cr_vars", clear
keep if bl_dm==0
gen antivegf=1
append using "$savedir\cataract\cr_vars"
replace antivegf = 0 if antivegf==.

dtable age_at_index i.gender i.eth5 i.bmi_cat i.smokstatus yrs_dm i.dm_type_f i.index_acarbose-index_nsaid i.bl_retinopathy i.bl_eye_dis i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_ncm i.bl_neph_syndrome i.bl_rrt i.bl_stroke i.bl_zoster i.bl_hf i.bl_pad i.bl_copd i.bl_neuropathy i.bl_amputation  prior_egfr i.prior_egfr_cat i.prior_ckd i.prior_acr_cat tot_appts_yr_prior scr_meas_yr_prior, by(antivegf) export($projdir\output\no_dm_table_1.xlsx, replace)

* Number of antivegf injections 
use "$savedir\antivegf\injections", clear
dtable number_antivegf_yr if bl_dm==0, by(yr_since_fu) export($projdir\output\no_dm_table_2.xlsx, replace)
dtable number_antivegf_yr if bl_dm==1, by(yr_since_fu) export($projdir\output\dm_table_2.xlsx, replace)

* Outcomes
* Diabetes groups
use "$savedir\antivegf\outcome_egfr_40", clear
merge 1:1 patid using "$savedir\antivegf\cr_vars", keepusing(bl_dm)
keep if bl_dm==1
gen antivegf=1
append using "$savedir\photocoag\outcome_egfr_40"
replace antivegf = 0 if antivegf==.

dtable i.egfr_40 i.egfr_40_sustained, by(antivegf)

stset end_egfr_40, id(patid) failure(egfr_40) enter(index_date)
strate antivegf, per(100000) 

stset end_egfr_40_sustained, id(patid) failure(egfr_40_sustained) enter(index_date)
strate antivegf, per(100000) 

* No diabetes groups
use "$savedir\antivegf\outcome_egfr_40", clear
merge 1:1 patid using "$savedir\antivegf\cr_vars", keepusing(bl_dm)
keep if bl_dm==0
gen antivegf=1
append using "$savedir\cataract\outcome_egfr_40"
replace antivegf = 0 if antivegf==.

dtable i.egfr_40 i.egfr_40_sustained, by(antivegf)

stset end_egfr_40, id(patid) failure(egfr_40) enter(index_date)
strate antivegf, per(100000) 

stset end_egfr_40_sustained, id(patid) failure(egfr_40_sustained) enter(index_date)
strate antivegf, per(100000) 

******************

