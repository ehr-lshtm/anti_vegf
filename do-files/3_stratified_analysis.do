/*==============================================================================
DO FILE NAME:			3_stratified_analysis.do
DATE: 					18/02/2025
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

* Diabetes group
use "$savedir\an_dm_main_analysis_ps", clear

* There are small numbers of events/people in some strata therefore will drop some group and collapse some categories
* Nephrotic syndrome too small numbers so not included in this analysis
* IMD - drop those with missing IMD 
* Ethnicity - combine non-White groups
* Albuminuria - often small numbers in A2 and A3 - combine
file open tablecontent using "$projdir/output/primary_cox_models_stratified.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Strata") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
* redefining variables 
* collapse high eGFR categories 
gen prior_egfr_cat_3 = prior_egfr_cat
replace prior_egfr_cat_3 = 30 if prior_egfr_cat<30
* Create binary ACR variable 
gen prior_acr_cat_bin = prior_acr_cat 
replace prior_acr_cat_bin = 1 if prior_acr_cat==2
* Ethnicity white vs non-white 
gen eth_bin = eth5 
replace eth_bin = 1 if eth5>1 & eth5<5
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	* modelling eGFR categories
	stcox antivegf if prior_egfr_cat_3==60, vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_`a'_egfr_normal", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_egfr_normal", replace) idstr("output_`a'_egfr_normal") 
	stcox antivegf if prior_egfr_cat_3==45, vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_`a'_egfr_mild", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_egfr_mild", replace) idstr("output_`a'_egfr_mild")
	stcox antivegf if prior_egfr_cat_3==30, vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_`a'_egfr_moderate", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_egfr_moderate", replace) idstr("output_`a'_egfr_moderate") 
	* Write to file - egfr normal
	count if antivegf==0 & prior_egfr_cat_3==60
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==60
	local event_0 = r(N)
	bysort antivegf prior_egfr_cat_3: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==60
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR normal") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==60
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==60
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==60
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'_egfr_normal" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* egfr mild
	count if antivegf==0 & prior_egfr_cat_3==45
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==45
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==45
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR mild") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==45
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==45
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==45
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR mild") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'_egfr_mild" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* egfr moderate
	count if antivegf==0 & prior_egfr_cat_3==30
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==30
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==30
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR moderate") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==30
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==30
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==30
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eGFR moderate") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'_egfr_moderate" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
	
	* modelling ACR categories
	stcox antivegf if prior_acr_cat_bin==0, vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_`a'_acr_normal", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_acr_normal", replace) idstr("output_`a'_acr_normal") 
	stcox antivegf if prior_acr_cat_bin==1, vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_`a'_acr_high", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_acr_high", replace) idstr("output_`a'_acr_high")
	* Write to file - acr normal
	count if antivegf==0 & prior_acr_cat_bin==0
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_acr_cat_bin==0
	local event_0 = r(N)
	bysort antivegf prior_acr_cat_bin: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0 & prior_acr_cat_bin==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("ACR normal") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_acr_cat_bin==0
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_acr_cat_bin==0
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_acr_cat_bin==0
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("ACR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'_acr_normal" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* acr high
	count if antivegf==0 & prior_acr_cat_bin==1
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_acr_cat_bin==1
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_acr_cat_bin==1
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("ACR high") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_acr_cat_bin==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_acr_cat_bin==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_acr_cat_bin==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("`a'") _tab ("ACR high") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_`a'_acr_high" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
	
	* IMD
	preserve 
	drop if imd==0 
	forvalues i=1/5 {
		stcox antivegf if imd==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_dm_`a'_imd_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_imd_`i'", replace) idstr("output_`a'_imd_`i'") 
		
		* Write to file 
		count if antivegf==0 & imd==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & imd==`i'
		local event_0 = r(N)
		bysort antivegf imd: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & imd==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("IMD `i'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & imd==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & imd==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & imd==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("IMD `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_dm_`a'_imd_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
	restore 
	
	* Ethnicity 
	foreach i in 0 1 5 {
		stcox antivegf if eth_bin==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_dm_`a'_eth_bin_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_eth_bin_`i'", replace) idstr("output_`a'_eth_bin_`i'") 
		
		* Write to file 
		count if antivegf==0 & eth_bin==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & eth_bin==`i'
		local event_0 = r(N)
		bysort antivegf eth_bin: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & eth_bin==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eth_bin `i'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & eth_bin==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & eth_bin==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & eth_bin==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("eth_bin `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_dm_`a'_eth_bin_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
		* Sex
	forvalues i=1/2 {
		stcox antivegf if gender==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_dm_`a'_gender_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_`a'_gender_`i'", replace) idstr("output_`a'_gender_`i'") 
		
		* Write to file 
		count if antivegf==0 & gender==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & gender==`i'
		local event_0 = r(N)
		bysort antivegf gender: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & gender==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("gender `i'") _tab ("Photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & gender==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & gender==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & gender==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("Diabetes") _tab ("`a'") _tab ("gender `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_dm_`a'_gender_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
}

* Hypertension outcome 
drop if bl_hypertension==1
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)

* modelling eGFR categories
stcox antivegf if prior_egfr_cat_3==60, vce(robust) 
estimates save "$projdir/output/tempdata/output_dm_event_hypertension_egfr_normal", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_egfr_normal", replace) idstr("output_event_hypertension_egfr_normal") 
stcox antivegf if prior_egfr_cat_3==45, vce(robust) 
estimates save "$projdir/output/tempdata/output_dm_event_hypertension_egfr_mild", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_egfr_mild", replace) idstr("output_event_hypertension_egfr_mild")
stcox antivegf if prior_egfr_cat_3==30, vce(robust) 
estimates save "$projdir/output/tempdata/output_dm_event_hypertension_egfr_moderate", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_egfr_moderate", replace) idstr("output_event_hypertension_egfr_moderate") 
* Write to file - egfr normal
count if antivegf==0 & prior_egfr_cat_3==60
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==60
local event_0 = r(N)
bysort antivegf prior_egfr_cat_3: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==60
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR normal") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==60
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==60
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==60
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_egfr_normal" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* egfr mild
count if antivegf==0 & prior_egfr_cat_3==45
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==45
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==45
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR mild") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==45
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==45
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==45
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR mild") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_egfr_mild" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* egfr moderate
count if antivegf==0 & prior_egfr_cat_3==30
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==30
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==30
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR moderate") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==30
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==30
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==30
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eGFR moderate") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_egfr_moderate" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* modelling ACR categories
stcox antivegf if prior_acr_cat_bin==0, vce(robust) 
estimates save "$projdir/output/tempdata/output_dm_event_hypertension_acr_normal", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_acr_normal", replace) idstr("output_event_hypertension_acr_normal") 
stcox antivegf if prior_acr_cat_bin==1, vce(robust) 
estimates save "$projdir/output/tempdata/output_dm_event_hypertension_acr_high", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_acr_high", replace) idstr("output_event_hypertension_acr_high")
* Write to file - acr normal
count if antivegf==0 & prior_acr_cat_bin==0
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_acr_cat_bin==0
local event_0 = r(N)
bysort antivegf prior_acr_cat_bin: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0 & prior_acr_cat_bin==0
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("ACR normal") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_acr_cat_bin==0
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_acr_cat_bin==0
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_acr_cat_bin==0
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("ACR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_acr_normal" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* acr high
count if antivegf==0 & prior_acr_cat_bin==1
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_acr_cat_bin==1
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_acr_cat_bin==1
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("ACR high") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_acr_cat_bin==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_acr_cat_bin==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_acr_cat_bin==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("ACR high") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_acr_high" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* IMD
preserve 
drop if imd==0 
forvalues i=1/5 {
	stcox antivegf if imd==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_event_hypertension_imd_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_imd_`i'", replace) idstr("output_event_hypertension_imd_`i'") 
	
	* Write to file 
	count if antivegf==0 & imd==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & imd==`i'
	local event_0 = r(N)
	bysort antivegf imd: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & imd==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("IMD `i'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & imd==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & imd==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & imd==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("IMD `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_imd_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}
restore 

* Ethnicity 
foreach i in 0 1 5 {
	stcox antivegf if eth_bin==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_event_hypertension_eth_bin_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_eth_bin_`i'", replace) idstr("output_event_hypertension_eth_bin_`i'") 
	
	* Write to file 
	count if antivegf==0 & eth_bin==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & eth_bin==`i'
	local event_0 = r(N)
	bysort antivegf eth_bin: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & eth_bin==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eth_bin `i'") _tab ("photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & eth_bin==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & eth_bin==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & eth_bin==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("eth_bin `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_eth_bin_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}
	* Sex
forvalues i=1/2 {
	stcox antivegf if gender==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_dm_event_hypertension_gender_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_dm_event_hypertension_gender_`i'", replace) idstr("output_event_hypertension_gender_`i'") 
	
	* Write to file 
	count if antivegf==0 & gender==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & gender==`i'
	local event_0 = r(N)
	bysort antivegf gender: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & gender==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("gender `i'") _tab ("Photocoagulation") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & gender==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & gender==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & gender==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("Diabetes") _tab ("event_hypertension") _tab ("gender `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_dm_event_hypertension_gender_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}


* No diabetes group*************************************************************************************************
use "$savedir\an_nodm_main_analysis_ps", clear

* redefining variables 
* collapse high eGFR categories 
gen prior_egfr_cat_3 = prior_egfr_cat
replace prior_egfr_cat_3 = 30 if prior_egfr_cat<30
* Create binary ACR variable 
gen prior_acr_cat_bin = prior_acr_cat 
replace prior_acr_cat_bin = 1 if prior_acr_cat==2
* Ethnicity white vs non-white 
gen eth_bin = eth5 
replace eth_bin = 1 if eth5>1 & eth5<5
* Fit weighted Cox regression w/ robust standard errors
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	* modelling eGFR categories
	stcox antivegf if prior_egfr_cat_3==60, vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_`a'_egfr_normal", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_egfr_normal", replace) idstr("output_`a'_egfr_normal") 
	stcox antivegf if prior_egfr_cat_3==45, vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_`a'_egfr_mild", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_egfr_mild", replace) idstr("output_`a'_egfr_mild")
	stcox antivegf if prior_egfr_cat_3==30, vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_`a'_egfr_moderate", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_egfr_moderate", replace) idstr("output_`a'_egfr_moderate") 
	* Write to file - egfr normal
	count if antivegf==0 & prior_egfr_cat_3==60
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==60
	local event_0 = r(N)
	bysort antivegf prior_egfr_cat_3: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==60
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR normal") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==60
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==60
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==60
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'_egfr_normal" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* egfr mild
	count if antivegf==0 & prior_egfr_cat_3==45
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==45
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==45
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR mild") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==45
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==45
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==45
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR mild") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'_egfr_mild" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* egfr moderate
	count if antivegf==0 & prior_egfr_cat_3==30
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_egfr_cat_3==30
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_egfr_cat_3==30
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR moderate") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_egfr_cat_3==30
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_egfr_cat_3==30
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_egfr_cat_3==30
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eGFR moderate") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'_egfr_moderate" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
	
	* modelling ACR categories
	stcox antivegf if prior_acr_cat_bin==0, vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_`a'_acr_normal", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_acr_normal", replace) idstr("output_`a'_acr_normal") 
	stcox antivegf if prior_acr_cat_bin==1, vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_`a'_acr_high", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_acr_high", replace) idstr("output_`a'_acr_high")
	* Write to file - acr normal
	count if antivegf==0 & prior_acr_cat_bin==0
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_acr_cat_bin==0
	local event_0 = r(N)
	bysort antivegf prior_acr_cat_bin: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0 & prior_acr_cat_bin==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("ACR normal") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_acr_cat_bin==0
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_acr_cat_bin==0
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_acr_cat_bin==0
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("ACR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'_acr_normal" 
    cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	* acr high
	count if antivegf==0 & prior_acr_cat_bin==1
	local denominator_0 = r(N)
	count if `a'==1 & antivegf==0 & prior_acr_cat_bin==1
	local event_0 = r(N)
	sum total_fu_`a' if antivegf==0 & prior_acr_cat_bin==1
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("ACR high") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & prior_acr_cat_bin==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1 & prior_acr_cat_bin==1
	local event_1 = r(N)
	sum total_fu_`a' if antivegf==1 & prior_acr_cat_bin==1
	local person_mth_1 = r(mean)/30
    local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("ACR high") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_`a'_acr_high" 
	cap lincom antivegf, eform
    file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
    cap estimates clear
	drop total_fu_`a'
	
	* IMD
	preserve 
	drop if imd==0 
	forvalues i=1/5 {
		stcox antivegf if imd==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_nodm_`a'_imd_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_imd_`i'", replace) idstr("output_`a'_imd_`i'") 
		
		* Write to file 
		count if antivegf==0 & imd==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & imd==`i'
		local event_0 = r(N)
		bysort antivegf imd: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & imd==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("IMD `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & imd==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & imd==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & imd==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("IMD `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_nodm_`a'_imd_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
	restore 
	
	* Ethnicity 
	foreach i in 0 1 5 {
		stcox antivegf if eth_bin==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_nodm_`a'_eth_bin_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_eth_bin_`i'", replace) idstr("output_`a'_eth_bin_`i'") 
		
		* Write to file 
		count if antivegf==0 & eth_bin==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & eth_bin==`i'
		local event_0 = r(N)
		bysort antivegf eth_bin: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & eth_bin==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eth_bin `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & eth_bin==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & eth_bin==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & eth_bin==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("eth_bin `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_nodm_`a'_eth_bin_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
	* Sex
	forvalues i=1/2 {
		stcox antivegf if gender==`i', vce(robust) 
		estimates save "$projdir/output/tempdata/output_nodm_`a'_gender_`i'", replace 
		eststo model1
		parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_`a'_gender_`i'", replace) idstr("output_`a'_gender_`i'") 
		
		* Write to file 
		count if antivegf==0 & gender==`i'
		local denominator_0 = r(N)
		count if `a'==1 & antivegf==0 & gender==`i'
		local event_0 = r(N)
		bysort antivegf gender: egen total_fu_`a' = total(_t)
		sum total_fu_`a' if antivegf==0 & gender==`i'
		local person_mth_0 = r(mean)/30
		local rate_0 = 100000*(`event_0'/`person_mth_0')
		di `rate_0'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("gender `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
		file write tablecontent ("1.00") _tab _tab ("1.00")  _n
		count if antivegf==1 & gender==`i'
		local denominator_1 = r(N)
		count if `a'==1 & antivegf==1 & gender==`i'
		local event_1 = r(N)
		sum total_fu_`a' if antivegf==1 & gender==`i'
		local person_mth_1 = r(mean)/30
		local rate_1 = 100000*(`event_1'/`person_mth_1')
		di `rate_1'
		file write tablecontent ("No diabetes") _tab ("`a'") _tab ("gender `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
		cap estimates use "$projdir/output/tempdata/output_nodm_`a'_gender_`i'" 
		cap lincom antivegf, eform
		file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
		drop total_fu_`a'
	}
	
}

* Hypertension outcome 
drop if bl_hypertension==1
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)

* modelling eGFR categories
stcox antivegf if prior_egfr_cat_3==60, vce(robust) 
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_normal", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_egfr_normal", replace) idstr("output_event_hypertension_egfr_normal") 
stcox antivegf if prior_egfr_cat_3==45, vce(robust) 
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_mild", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_egfr_mild", replace) idstr("output_event_hypertension_egfr_mild")
stcox antivegf if prior_egfr_cat_3==30, vce(robust) 
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_moderate", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_egfr_moderate", replace) idstr("output_event_hypertension_egfr_moderate") 
* Write to file - egfr normal
count if antivegf==0 & prior_egfr_cat_3==60
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==60
local event_0 = r(N)
bysort antivegf prior_egfr_cat_3: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==60
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR normal") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==60
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==60
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==60
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_normal" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* egfr mild
count if antivegf==0 & prior_egfr_cat_3==45
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==45
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==45
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR mild") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==45
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==45
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==45
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR mild") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_mild" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* egfr moderate
count if antivegf==0 & prior_egfr_cat_3==30
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_egfr_cat_3==30
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_egfr_cat_3==30
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR moderate") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_egfr_cat_3==30
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_egfr_cat_3==30
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_egfr_cat_3==30
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eGFR moderate") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_egfr_moderate" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* modelling ACR categories
stcox antivegf if prior_acr_cat_bin==0, vce(robust) 
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_acr_normal", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_acr_normal", replace) idstr("output_event_hypertension_acr_normal") 
stcox antivegf if prior_acr_cat_bin==1, vce(robust) 
estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_acr_high", replace 
eststo model1
parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_acr_high", replace) idstr("output_event_hypertension_acr_high")
* Write to file - acr normal
count if antivegf==0 & prior_acr_cat_bin==0
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_acr_cat_bin==0
local event_0 = r(N)
bysort antivegf prior_acr_cat_bin: egen total_fu_event_hypertension = total(_t)
sum total_fu_event_hypertension if antivegf==0 & prior_acr_cat_bin==0
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("ACR normal") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_acr_cat_bin==0
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_acr_cat_bin==0
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_acr_cat_bin==0
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("ACR normal") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_acr_normal" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
* acr high
count if antivegf==0 & prior_acr_cat_bin==1
local denominator_0 = r(N)
count if event_hypertension==1 & antivegf==0 & prior_acr_cat_bin==1
local event_0 = r(N)
sum total_fu_event_hypertension if antivegf==0 & prior_acr_cat_bin==1
local person_mth_0 = r(mean)/30
local rate_0 = 100000*(`event_0'/`person_mth_0')
di `rate_0'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("ACR high") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
file write tablecontent ("1.00") _tab _tab ("1.00")  _n
count if antivegf==1 & prior_acr_cat_bin==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1 & prior_acr_cat_bin==1
local event_1 = r(N)
sum total_fu_event_hypertension if antivegf==1 & prior_acr_cat_bin==1
local person_mth_1 = r(mean)/30
local rate_1 = 100000*(`event_1'/`person_mth_1')
di `rate_1'
file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("ACR high") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_acr_high" 
cap lincom antivegf, eform
file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
cap estimates clear
drop total_fu_event_hypertension

* IMD
preserve 
drop if imd==0 
forvalues i=1/5 {
	stcox antivegf if imd==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_imd_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_imd_`i'", replace) idstr("output_event_hypertension_imd_`i'") 
	
	* Write to file 
	count if antivegf==0 & imd==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & imd==`i'
	local event_0 = r(N)
	bysort antivegf imd: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & imd==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("IMD `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & imd==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & imd==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & imd==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("IMD `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_imd_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}
restore 

* Ethnicity 
foreach i in 0 1 5 {
	stcox antivegf if eth_bin==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_eth_bin_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_eth_bin_`i'", replace) idstr("output_event_hypertension_eth_bin_`i'") 
	
	* Write to file 
	count if antivegf==0 & eth_bin==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & eth_bin==`i'
	local event_0 = r(N)
	bysort antivegf eth_bin: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & eth_bin==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eth_bin `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & eth_bin==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & eth_bin==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & eth_bin==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("eth_bin `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_eth_bin_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}
* Sex
forvalues i=1/2 {
	stcox antivegf if gender==`i', vce(robust) 
	estimates save "$projdir/output/tempdata/output_nodm_event_hypertension_gender_`i'", replace 
	eststo model1
	parmest, label eform format(estimate p lb ub) saving("$projdir/output/tempdata/output_nodm_event_hypertension_gender_`i'", replace) idstr("output_event_hypertension_gender_`i'") 
	
	* Write to file 
	count if antivegf==0 & gender==`i'
	local denominator_0 = r(N)
	count if event_hypertension==1 & antivegf==0 & gender==`i'
	local event_0 = r(N)
	bysort antivegf gender: egen total_fu_event_hypertension = total(_t)
	sum total_fu_event_hypertension if antivegf==0 & gender==`i'
	local person_mth_0 = r(mean)/30
	local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("gender `i'") _tab ("Cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
	file write tablecontent ("1.00") _tab _tab ("1.00")  _n
	count if antivegf==1 & gender==`i'
	local denominator_1 = r(N)
	count if event_hypertension==1 & antivegf==1 & gender==`i'
	local event_1 = r(N)
	sum total_fu_event_hypertension if antivegf==1 & gender==`i'
	local person_mth_1 = r(mean)/30
	local rate_1 = 100000*(`event_1'/`person_mth_1')
	di `rate_1'
	file write tablecontent ("No diabetes") _tab ("event_hypertension") _tab ("gender `i'") _tab ("antivegf") _tab (`denominator_1') _tab (`event_1') _tab %10.0f (`person_mth_1') _tab %3.2f (`rate_1') _tab 
	cap estimates use "$projdir/output/tempdata/output_nodm_event_hypertension_gender_`i'" 
	cap lincom antivegf, eform
	file write tablecontent  %4.2f (r(estimate)) _tab ("(") %4.2f (r(lb)) (" - ") %4.2f (r(ub)) (")") _tab %4.2f (r(lb)) _tab %4.2f (r(ub)) _tab _n
	drop total_fu_event_hypertension
}

file close tablecontent