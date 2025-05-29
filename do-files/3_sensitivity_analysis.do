/*==============================================================================
DO FILE NAME:			3_sensitivity_analysis.do
DATE: 					19/02/2025
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	runs main analysis

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\3_stratified_analysis_$S_DATE.log", append

** Sensitivity analysis 1: Per protocol analysis
* Explore starting other treatment 
use "$savedir\antivegf\cr_study_pop", clear
merge 1:1 patid using "$savedir\first_code_inclusions", keep(match)
egen total_procedures = rownonmiss( first_code_antivegf first_code_photocoag first_code_cataract)
tab total_procedures
gen time_cataract = first_code_cataract - first_code_antivegf 
gen time_photocoag = first_code_photocoag - first_code_antivegf 
bys bl_dm: sum time*, d
use "$savedir\photocoag\cr_study_pop", clear
merge 1:1 patid using "$savedir\first_code_inclusions", keep(match)
egen total_procedures = rownonmiss( first_code_antivegf first_code_photocoag first_code_cataract)
tab total_procedures
gen time_cataract = first_code_cataract - first_code_photocoag
gen time_antivegf = first_code_antivegf - first_code_photocoag
sum time*
use "$savedir\cataract\cr_study_pop", clear
merge 1:1 patid using "$savedir\first_code_inclusions", keep(match)
egen total_procedures = rownonmiss( first_code_antivegf first_code_photocoag first_code_cataract)
tab total_procedures
gen time_photocoag = first_code_photocoag - first_code_cataract
gen time_antivegf = first_code_antivegf - first_code_cataract
sum time*

use "$savedir\antivegf\injections", clear
keep patid last_injection gap_6_injections first_gap_date 
duplicates drop 
count 
codebook patid 
tempfile tempfile 
save `tempfile'

* Diabetes group
use "$savedir\an_dm_main_analysis_ps", clear
merge 1:1 patid using `tempfile'
drop if _merge==2 
tab _merge antivegf 
drop if _merge==1 & antivegf==1
drop _merge 
* censor 6 months after last injection or at first 6 month gap 
gen last_injection_6 = first_gap_date + 183
gen last_ever_injection_6 = last_injection + 183 
count if last_injection_6<end_fu
count if last_ever_injection_6<end_fu & last_ever_injection_6 < last_injection_6
* Create new end of follow-up where exposure ends at end of antivegf exposure 
gen end_fu_pp = end_fu 
replace end_fu_pp = last_injection_6 if last_injection_6<end_fu 
replace end_fu_pp = last_ever_injection_6 if last_ever_injection_6<end_fu_pp

* Update end of follow-up for unexposed to end if start antivegf 
merge 1:1 patid using "$savedir\first_code_inclusions", keep(match) keepusing(first_code_antivegf)
replace end_fu_pp = first_code_antivegf if antivegf==0 & first_code_antivegf<end_fu
sum end_fu*, d

local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	gen `b'_pp = end_fu_pp < `b'
	tab `b'_pp `a'
	gen time_`a'_pp = `b' - end_fu_pp
	bys bl_dm: sum time_`a'_pp if `b'_pp==1, d
	replace `b' = end_fu_pp if `b'_pp==1 
	replace `a' = 0 if `b'_pp==1 & `a'==1
}

dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_cvd_death i.event_zoster, by(antivegf)

* If imbalnaced - linear variables - quadratic/splines, or interaction terms

file open tablecontent using "$projdir/output/primary_cox_models_sensitivity_1.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("crude hr") _tab ("crude ci") _tab ("crude lci") _tab ("crude uci") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b', failure(`a') origin(index_date) enter(index_date) id(patid)
	sum _t, d 
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_dm_crude_`a'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_crude_`a'", replace) idstr("output_crude_`a'") 
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_dm_`a'", replace 
	eststo model2
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'", replace) idstr("output_`a'") 
	count if antivegf==0
	local denominator_0 = r(N)
	count if _d==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if _d==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_crude_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab
    cap estimates clear
	cap estimates use "$projdir/output/tempdata/output_dm_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
}

* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
drop if bl_hypertension==1
stset end_hypertension, failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_dm_crude_event_hypertension", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_crude_event_hypertension", replace) idstr("output_crude_event_hypertension") 
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_dm_event_hypertension", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension", replace) idstr("output_event_hypertension") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if _d==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if _d==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_crude_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab
cap estimates clear
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* No diabetes

use "$savedir\an_nodm_main_analysis_ps", clear
merge 1:1 patid using `tempfile'
drop if _merge==2 
tab _merge antivegf 
drop if _merge==1 & antivegf==1
drop _merge 
* censor 6 months after last injection or at first 6 month gap 
gen last_injection_6 = first_gap_date + 183
gen last_ever_injection_6 = last_injection + 183 
count if last_injection_6<end_fu
count if last_ever_injection_6<end_fu & last_ever_injection_6 < last_injection_6
* Create new end of follow-up where exposure ends at end of antivegf exposure 
gen end_fu_pp = end_fu 
replace end_fu_pp = last_injection_6 if last_injection_6<end_fu 
replace end_fu_pp = last_ever_injection_6 if last_ever_injection_6<end_fu_pp

* Update end of follow-up for unexposed to end if start antivegf 
merge 1:1 patid using "$savedir\first_code_inclusions", keep(match) keepusing(first_code_antivegf)
replace end_fu_pp = first_code_antivegf if antivegf==0 & first_code_antivegf<end_fu
sum end_fu*, d

local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	gen `b'_pp = end_fu_pp < `b'
	replace `b' = end_fu_pp if `b'_pp==1 
	replace `a' = 0 if `b'_pp==1
}

dtable i.egfr_40 i.egfr_40_sustained i.acr_increased i.event_af i.event_hypertension i.event_kidney_failure i.event_mi i.event_neph_syndrome i.event_stroke i.event_hf i.event_pad i.event_zoster, by(antivegf)

* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b', failure(`a') origin(index_date) enter(index_date) id(patid)
	sum _t, d 
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_nodm_crude_`a'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_crude_`a'", replace) idstr("output_crude_`a'") 
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_nodm_`a'", replace 
	eststo model2
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'", replace) idstr("output_`a'") 
	count if antivegf==0
	local denominator_0 = r(N)
	count if _d==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if _d==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent  ("No diabetes") _tab ("`a'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_crude_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab
    cap estimates clear
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
}

* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
drop if bl_hypertension==1
stset end_hypertension, failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_nodm_crude_event_hypertension", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_crude_event_hypertension", replace) idstr("output_crude_event_hypertension") 
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension", replace 
eststo model2
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension", replace) idstr("output_event_hypertension") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if _d==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if _d==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent  ("No diabetes") _tab ("event_hypertension") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_crude_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab
cap estimates clear
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

file close tablecontent 

** Sensitivity analysis 2: Outcome definition 
* Determine median using measurements 3 months prior to index 
foreach grp in photocoag cataract {
	use "$savedir\\`grp'\eGFR", clear
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match)
	codebook patid 
	gen meas_time_index = obsdate - index_date
	gen meas_3_mth = meas_time_index>=-92 & meas_time_index<=0
	bys patid: egen tot_meas_3 = total(meas_3_mth)
	bys patid: egen median_3 = median(egfr) if meas_3_mth==1
	keep patid median_3
	duplicates drop 
	sum median_3, d
}
use "$savedir\antivegf\eGFR", clear
merge m:1 patid using "$savedir\antivegf\cr_study_pop", keepusing(index_date bl_dm) keep(match)
codebook patid 
gen meas_time_index = obsdate - index_date
gen meas_3_mth = meas_time_index>=-92 & meas_time_index<=0
bys patid: egen tot_meas_3 = total(meas_3_mth)
bys patid: egen median_3 = median(egfr) if meas_3_mth==1
keep patid median_3 bl_dm
duplicates drop 
bys bl_dm: sum median_3, d

** Sensitivity analysis 3: People with coded diabetic retinopathy only 
* Diabetes group
use "$savedir\an_dm_main_analysis_ps", clear
keep if bl_retinopathy==1

file open tablecontent using "$projdir/output/primary_cox_models_sensitivity_3.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_dm_`a'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'", replace) idstr("output_`a'") 
	count if antivegf==0
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
}

* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
drop if bl_hypertension==1
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_dm_event_hypertension", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension", replace) idstr("output_event_hypertension") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if event_hypertension==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

file close tablecontent

** Sensitivity analysis 4: Coded dialysis or kidney transplant only for kidney failure

* Diabetes group 
use "$savedir\an_dm_main_analysis_ps", clear
stset end_kidney_failure [pweight=att_weight], failure(event_kidney_failure) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_dm_event_kidney_failure", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_kidney_failure", replace) idstr("output_event_kidney_failure") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if event_kidney_failure==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_kidney_failure = total(_t)
sum total_fu_event_kidney_failure if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file open tablecontent using "$projdir/output/primary_cox_models_sensitivity_4.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
file write tablecontent ("Diabetes") _tab ("event_kidney_failure") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_kidney_failure==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_kidney_failure if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent  ("Diabetes") _tab ("event_kidney_failure") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_kidney_failure" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_kidney_failure

use "$savedir\an_nodm_main_analysis_ps", clear
stset end_kidney_failure [pweight=att_weight], failure(event_kidney_failure) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_nodm_event_kidney_failure", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_kidney_failure", replace) idstr("output_event_kidney_failure") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if event_kidney_failure==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_kidney_failure = total(_t)
sum total_fu_event_kidney_failure if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_kidney_failure") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_kidney_failure==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_kidney_failure if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent  ("No diabetes") _tab ("event_kidney_failure") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_kidney_failure" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_kidney_failure

file close tablecontent

* Sensitivity analysis 5: Stratifying on time: restricting follow-up to 0-6-months, 0-12 months, 0-24 months, 0-36 months, 0-48 months & 0-60 months.

use "$savedir\an_dm_main_analysis_ps", clear

file open tablecontent using "$projdir/output/primary_cox_models_sensitivity_5.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Time") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
* Fit weighted Cox regression w/ robust standard errors
foreach value in 183 365.25 730.5 1095.75 1461 1826.25 {
	local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
	local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
	local n: word count `event'
	forvalues i=1/`n' {
		local a : word `i' of `event'
		local b : word `i' of `end_date'
		stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid) exit(time index_date + `value')
		stcox antivegf, vce(robust)
		estimates save "$projdir/output/tempdata/output_dm_`a'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'", replace) idstr("output_`a'") 
		count if antivegf==0
		local denominator_0 = r(N)
		count if _d==1 & antivegf==0
		local event_0 = r(N)
		bysort antivegf: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("`value'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1
		local denominator_1 = r(N)
		count if _d==1 & antivegf==1
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("`value'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_dm_`a'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		cap estimates clear
		drop total_fu_`a'
	}

	* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
	drop if bl_hypertension==1
	stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid) exit(time index_date + `value')
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_dm_event_hypertension", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension", replace) idstr("output_event_hypertension") 
	count if antivegf==0
	local denominator_0 = r(N)
	di `denominator_0'
	count if _d==1 & antivegf==0
	local event_0 = r(N)
	di `event_0'
	bysort antivegf: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0
	local person_mth_0 = r(mean)/30
	di `person_mth_0'
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("`value'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if _d==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("`value'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	cap estimates clear
	drop total_fu_event_hypertension
}

* No diabetes

use "$savedir\an_nodm_main_analysis_ps", clear
foreach value in 183 365.25 730.5 1095.75 1461 1826.25 {
* Fit weighted Cox regression w/ robust standard errors
	local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
	local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
	local n: word count `event'
	forvalues i=1/`n' {
		local a : word `i' of `event'
		local b : word `i' of `end_date'
		stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid) exit(time index_date + `value')
		stcox antivegf, vce(robust)
		estimates save "$projdir/output/tempdata/output_nodm_`a'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'", replace) idstr("output_`a'") 
		count if antivegf==0
		local denominator_0 = r(N)
		count if _d==1 & antivegf==0
		local event_0 = r(N)
		bysort antivegf: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("`value'") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1
		local denominator_1 = r(N)
		count if _d==1 & antivegf==1
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent  ("No diabetes") _tab ("`a'") _tab ("`value'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_nodm_`a'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		cap estimates clear
		drop total_fu_`a'
	}

	* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
	drop if bl_hypertension==1
	stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)  exit(time index_date + `value')
		stcox antivegf, vce(robust)
		estimates save "$projdir/output/tempdata/output_nodm_event_hypertension", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension", replace) idstr("output_event_hypertension") 
	count if antivegf==0
	local denominator_0 = r(N)
	di `denominator_0'
	count if event_hypertension==1 & antivegf==0
	local event_0 = r(N)
	di `event_0'
	bysort antivegf: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0
	local person_mth_0 = r(mean)/30
	di `person_mth_0'
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("`value'") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent  ("No diabetes") _tab ("event_hypertension") _tab ("`value'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	cap estimates clear
	drop total_fu_event_hypertension
}

file close tablecontent

** Sensitivity analysis 6: People without missing ethnicity
* Diabetes group
use "$savedir\an_dm_main_analysis_ps", clear
keep if eth5!=5

file open tablecontent using "$projdir/output/primary_cox_models_sensitivity_6.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_dm_`a'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'", replace) idstr("output_`a'") 
	count if antivegf==0
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
}

* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
drop if bl_hypertension==1
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
stcox antivegf, vce(robust)
estimates save "$projdir/output/tempdata/output_dm_event_hypertension", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension", replace) idstr("output_event_hypertension") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if event_hypertension==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* No diabetes 
use "$savedir\an_nodm_main_analysis_ps", clear
keep if eth5!=5
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_nodm_`a'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'", replace) idstr("output_`a'") 
	count if antivegf==0
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent  ("No diabetes") _tab ("`a'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
}

* Hypertension outcome: Use only people without hypertension at index as can't have new diagnosis of hypertension 
drop if bl_hypertension==1
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)
	stcox antivegf, vce(robust)
	estimates save "$projdir/output/tempdata/output_nodm_event_hypertension", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension", replace) idstr("output_event_hypertension") 
count if antivegf==0
local denominator_0 = r(N)
di `denominator_0'
count if event_hypertension==1 & antivegf==0
local event_0 = r(N)
di `event_0'
bysort antivegf: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0
local person_mth_0 = r(mean)/30
di `person_mth_0'
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00") _tab _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent  ("No diabetes") _tab ("event_hypertension") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

file close tablecontent