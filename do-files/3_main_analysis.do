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

* Diabetes group
use "$savedir\an_dm_main_analysis", clear
* Model treatment allocation on the set of confounding variables  - want things that are related to outcome, not only related to exposure
logistic antivegf age_at_index i.gender i.eth5 ib1.imd i.smokstatus i.bmi_cat yrs_dm dm_type i.drug_dm_count yrs_retinopathy i.bl_amputation i.bl_neuropathy i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_neph_syndrome i.bl_stroke i.bl_hf i.bl_pad i.bl_copd i.bl_cancer yr_index i.index_statin i.index_acei i.index_antiplatelet i.index_arb i.index_betablocker i.index_ccb i.index_loop_diuretic i.index_mra i.index_oac i.index_otherantihypertensive i.index_ppi i.index_nsaid scr_meas_yr_prior tot_appts_yr_prior

* Estimate propensity scores
predict propensity

graph tw kdensity propensity if antivegf == 0 || kdensity propensity if antivegf == 1, leg(lab(1 "Photocoagulation") lab(2 "Antivegf"))
graph export "$projdir\output\ps_dm_all.jpg", replace
* Summarise distributions - 
bys antivegf: sum propensity, d
* Balance diagnostics
gen wts = 1/propensity if antivegf == 1
replace wts = 1/(1-propensity) if antivegf == 0

* Trim propensity score
egen min_ps_1 = min(propensity) if antivegf==1
egen min_ps_0 = min(propensity) if antivegf==0
egen max_ps_1 = max(propensity) if antivegf==1
egen max_ps_0 = max(propensity) if antivegf==0
foreach var in min_ps_1 min_ps_0 max_ps_1 max_ps_0 {
	sort `var'
	replace `var' = `var'[_n-1] if `var'==.
}

gen trim_min = max(min_ps_1, min_ps_0)
gen trim_max = min(max_ps_1, max_ps_0)

count if propensity<trim_min 
count if propensity>trim_max
* drop rather than generate new variable 
drop if (propensity<trim_min | propensity>trim_max)

* Calculate ATT weights 
gen att_weight = antivegf + (1-antivegf)*(propensity/(1-propensity))

* Without weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat yrs_dm dm_type drug_dm_count yrs_retinopathy bl_amputation bl_neuropathy bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs
* With weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat yrs_dm dm_type drug_dm_count yrs_retinopathy bl_amputation bl_neuropathy bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs wt(att_weight)

save "$savedir\an_dm_main_analysis_ps", replace
* If imbalnaced - linear variables - quadratic/splines, or interaction terms

file open tablecontent using "$projdir/output/primary_cox_models.txt", write text replace
file write tablecontent ("Group") _tab ("Outcome") _tab ("Exposure group") _tab ("denominator") _tab ("events") _tab ("total_person_mth") _tab ("Rate") _tab ("Crude hr") _tab ("Crude ci") _tab ("Crude lci") _tab ("Crude uci") _tab ("hr") _tab ("ci") _tab ("lci") _tab ("uci") _tab  _n
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
    file write tablecontent ("1.00") _tab _tab ("1.00") _tab ("1.00") _n
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

* Generate cumulative incidence curves
file open tablech using "$projdir/output/cumulative_hazard.txt", write text replace
	file write tablech ("Group")  ("Outcome") _tab ("Exposure") _tab ("Cumulative incidence") _tab ("95% confidence interval") _n
* Cumulative hazard plots 
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure end_mi end_neph_syndrome end_stroke end_hf end_pad end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	
	* Add cumulative incidence plots 
	* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
	* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
	stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
	summ _t
	local tmax=r(max)
	local tmaxplus1=r(max)+1

	range days 0 `tmax' `tmaxplus1'
	stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

	gen date = index_date + days
	format date %tddd_Month

	for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X
	
	*cumulative outcomes at last day of follow-up - write to file 
	
	file write tablech ("Diabetes") _tab ("`a'") _tab ("Photoocagulation") _tab 
	* cumulative outcome - photocoagulation 
	sum _at1 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 
	* cumulative outcome - Antivegf 
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
	sum _at2 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 

	*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
	********************* Need to sort out graph axes
	twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
					(rarea _at2_lci _at2_uci days, color(blue%25)) ///
					(line _at1 days, sort lcolor(red)) ///
					(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
					, legend(order(1 "Photocoagulation" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
					title("Time to `a'", justification(left) size(med) )  	   ///
					yscale(range(0, 1)) 											///
					ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
					xlabel(0 (500) 2700, labsize(small))				   				///			
					ytitle("Cumulative outcomes (%)", size(medsmall)) ///
					xtitle("days since index date", size(medsmall))      		///
					graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

	graph export "$projdir/output/adjcurv_dm_`a'.svg", as(svg) replace

	* Close window 
	graph close
	drop days-date
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
eststo model2
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
file write tablecontent ("1.00") _tab _tab ("1.00") _tab ("1.00") _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1
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


* Cumulative hazard plots 

stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)

* Add cumulative incidence plots 
* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
summ _t
local tmax=r(max)
local tmaxplus1=r(max)+1

range days 0 `tmax' `tmaxplus1'
stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

gen date = index_date + days
format date %tddd_Month

for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X

*cumulative outcomes at last day of follow-up - write to file 

file write tablech ("Diabetes") _tab ("event_hypertension") _tab ("Photoocagulation") _tab 
* cumulative outcome - photocoagulation 
sum _at1 if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at1_lci if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at1_uci if days==`tmax'
file write tablech (r(mean)) _tab _n 
* cumulative outcome - Antivegf 
file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
sum _at2 if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at2_lci if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at2_uci if days==`tmax'
file write tablech (r(mean)) _tab _n 

*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
********************* Need to sort out graph axes
twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
				(rarea _at2_lci _at2_uci days, color(blue%25)) ///
				(line _at1 days, sort lcolor(red)) ///
				(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
				, legend(order(1 "Photocoagulation" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
				title("Time to hypertension", justification(left) size(med) )  	   ///
				yscale(range(0, 1)) 											///
				ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
				xlabel(0 (500) 2700, labsize(small))				   				///			
				ytitle("Cumulative outcomes (%)", size(medsmall)) ///
				xtitle("days since index date", size(medsmall))      		///
				graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

graph export "$projdir/output/adjcurv_dm_event_hypertension.svg", as(svg) replace

* Close window 
graph close
drop days-date


/*use "$projdir/output/tempdata/output_dm_egfr_40", clear
foreach event in egfr_40_sustained acr_increased event_af event_hypertension event_kidney_failure event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster {
	append using "$projdir/output/tempdata/output_dm_`event'"
}
save "$projdir/output/cox_output_dm_primary", replace
export delimited using "$projdir/output/cox_output_dm_primary.csv", replace */




* No diabetes
use "$savedir\an_nodm_main_analysis", clear
* Model treatment allocation on the set of confounding variables  - want things that are related to outcome, not only related to exposure
logistic antivegf age_at_index i.gender i.eth5 i.imd i.smokstatus i.bmi_cat i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_neph_syndrome i.bl_stroke i.bl_hf i.bl_pad i.bl_copd i.bl_cancer yr_index i.index_statin i.index_acei i.index_antiplatelet i.index_arb i.index_betablocker i.index_ccb i.index_loop_diuretic i.index_mra i.index_oac i.index_otherantihypertensive i.index_ppi i.index_nsaid scr_meas_yr_prior tot_appts_yr_prior

* Estimate propensity scores
predict propensity

graph tw kdensity propensity if antivegf == 0 || kdensity propensity if antivegf == 1, leg(lab(1 "Cataract") lab(2 "Antivegf"))
graph export "$projdir\output\ps_nodm_all.jpg", replace
* Summarise distributions - 
bys antivegf: sum propensity, d
* Balance diagnostics
gen wts = 1/propensity if antivegf == 1
replace wts = 1/(1-propensity) if antivegf == 0

* Trim propensity score
egen min_ps_1 = min(propensity) if antivegf==1
egen min_ps_0 = min(propensity) if antivegf==0
egen max_ps_1 = max(propensity) if antivegf==1
egen max_ps_0 = max(propensity) if antivegf==0
foreach var in min_ps_1 min_ps_0 max_ps_1 max_ps_0 {
	sort `var'
	replace `var' = `var'[_n-1] if `var'==.
}

gen trim_min = max(min_ps_1, min_ps_0)
gen trim_max = min(max_ps_1, max_ps_0)

count if propensity<trim_min 
count if propensity>trim_max

* drop rather than generate new variable 
drop if (propensity<trim_min | propensity>trim_max)

* Calculate ATT weights 
gen att_weight = antivegf + (1-antivegf)*(propensity/(1-propensity))

* Without weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs
* With weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs wt(att_weight)

save "$savedir\an_nodm_main_analysis_ps", replace
* If imbalnaced - linear variables - quadratic/splines, or interaction terms

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
	count if `a'==1 & antivegf==0
	local event_0 = r(N)
	bysort antivegf: egen total_fu_`a' = total(_t)
	sum total_fu_`a' if antivegf==0
	local person_mth_0 = r(mean)/30
    local rate_0 = 100000*(`event_0'/`person_mth_0')
	di `rate_0'
	file write tablecontent ("No diabetes") _tab ("`a'") _tab ("cataract") _tab (`denominator_0') _tab (`event_0') _tab %10.0f (`person_mth_0') _tab %3.2f (`rate_0') _tab 
    file write tablecontent ("1.00") _tab _tab ("1.00") _tab ("1.00") _n
	count if antivegf==1
	local denominator_1 = r(N)
	count if `a'==1 & antivegf==1
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

* Cumulative hazard plots 
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure end_mi end_neph_syndrome end_stroke end_hf end_pad end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	
	* Add cumulative incidence plots 
	* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
	* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
	stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
	summ _t
	local tmax=r(max)
	local tmaxplus1=r(max)+1

	range days 0 `tmax' `tmaxplus1'
	stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

	gen date = index_date + days
	format date %tddd_Month

	for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X
	
	*cumulative outcomes at last day of follow-up - write to file 
	
	file write tablech ("No diabetes") _tab ("`a'") _tab ("Cataract") _tab 
	* cumulative outcome - photocoagulation 
	sum _at1 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 
	* cumulative outcome - Antivegf 
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
	sum _at2 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 

	*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
	********************* Need to sort out graph axes
	twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
					(rarea _at2_lci _at2_uci days, color(blue%25)) ///
					(line _at1 days, sort lcolor(red)) ///
					(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
					, legend(order(1 "Cataract" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
					title("Time to `a'", justification(left) size(med) )  	   ///
					yscale(range(0, 1)) 											///
					ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
					xlabel(0 (500) 2700, labsize(small))				   				///			
					ytitle("Cumulative outcomes (%)", size(medsmall)) ///
					xtitle("days since index date", size(medsmall))      		///
					graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

	graph export "$projdir/output/adjcurv_nodm_`a'.svg", as(svg) replace

	* Close window 
	graph close
	drop days-date
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
file write tablecontent ("1.00") _tab _tab ("1.00") _tab ("1.00") _tab _tab _tab ("1.00") _n
count if antivegf==1
local denominator_1 = r(N)
count if event_hypertension==1 & antivegf==1
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

* Cumulative hazard plots 
stset end_hypertension [pweight=att_weight], failure(event_hypertension) origin(index_date) enter(index_date) id(patid)

* Add cumulative incidence plots 
* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
summ _t
local tmax=r(max)
local tmaxplus1=r(max)+1

range days 0 `tmax' `tmaxplus1'
stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

gen date = index_date + days
format date %tddd_Month

for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X

*cumulative outcomes at last day of follow-up - write to file 

file write tablech ("No diabetes") _tab ("event_hypertension") _tab ("Cataract") _tab 
* cumulative outcome - photocoagulation 
sum _at1 if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at1_lci if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at1_uci if days==`tmax'
file write tablech (r(mean)) _tab _n 
* cumulative outcome - Antivegf 
file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
sum _at2 if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at2_lci if days==`tmax'
file write tablech (r(mean)) _tab 
sum _at2_uci if days==`tmax'
file write tablech (r(mean)) _tab _n 

*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
********************* Need to sort out graph axes
twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
				(rarea _at2_lci _at2_uci days, color(blue%25)) ///
				(line _at1 days, sort lcolor(red)) ///
				(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
				, legend(order(1 "Cataract" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
				title("Time to hypertension", justification(left) size(med) )  	   ///
				yscale(range(0, 1)) 											///
				ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
				xlabel(0 (500) 2700, labsize(small))				   				///			
				ytitle("Cumulative outcomes (%)", size(medsmall)) ///
				xtitle("days since index date", size(medsmall))      		///
				graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

graph export "$projdir/output/adjcurv_nodm_event_hypertension.svg", as(svg) replace

* Close window 
graph close
drop days-date

file close tablech

/* cumulative incidence plots
 
use "$savedir\an_dm_main_analysis", clear
logistic antivegf age_at_index i.gender i.eth5 ib1.imd i.smokstatus i.bmi_cat yrs_dm dm_type i.drug_dm_count yrs_retinopathy i.bl_amputation i.bl_neuropathy i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_neph_syndrome i.bl_stroke i.bl_hf i.bl_pad i.bl_copd i.bl_cancer yr_index i.index_statin i.index_acei i.index_antiplatelet i.index_arb i.index_betablocker i.index_ccb i.index_loop_diuretic i.index_mra i.index_oac i.index_otherantihypertensive i.index_ppi i.index_nsaid scr_meas_yr_prior tot_appts_yr_prior

* Estimate propensity scores
predict propensity

graph tw kdensity propensity if antivegf == 0 || kdensity propensity if antivegf == 1
*graph export "$projdir\output\ps_dm_all.jpg", replace
* Summarise distributions - 
bys antivegf: sum propensity, d
* Balance diagnostics
gen wts = 1/propensity if antivegf == 1
replace wts = 1/(1-propensity) if antivegf == 0

* Trim propensity score
egen min_ps_1 = min(propensity) if antivegf==1
egen min_ps_0 = min(propensity) if antivegf==0
egen max_ps_1 = max(propensity) if antivegf==1
egen max_ps_0 = max(propensity) if antivegf==0
foreach var in min_ps_1 min_ps_0 max_ps_1 max_ps_0 {
	sort `var'
	replace `var' = `var'[_n-1] if `var'==.
}

gen trim_min = max(min_ps_1, min_ps_0)
gen trim_max = min(max_ps_1, max_ps_0)

count if propensity<trim_min 
count if propensity>trim_max
* drop rather than generate new variable 
drop if (propensity<trim_min | propensity>trim_max)

* Calculate ATT weights 
gen att_weight = antivegf + (1-antivegf)*(propensity/(1-propensity))

* Without weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat yrs_dm dm_type drug_dm_count yrs_retinopathy bl_amputation bl_neuropathy bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs
* With weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat yrs_dm dm_type drug_dm_count yrs_retinopathy bl_amputation bl_neuropathy bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs wt(att_weight)

file open tablech using "$projdir/output/cumulative_hazard.txt", write text replace
	file write tablech ("Group")  ("Outcome") _tab ("Exposure") _tab ("Cumulative incidence") _tab ("95% confidence interval") _n
* Cumulative hazard plots 
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure end_mi end_neph_syndrome end_stroke end_hf end_pad end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	
	* Add cumulative incidence plots 
	* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
	* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
	stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
	summ _t
	local tmax=r(max)
	local tmaxplus1=r(max)+1

	range days 0 `tmax' `tmaxplus1'
	stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

	gen date = index_date + days
	format date %tddd_Month

	for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X
	
	*cumulative outcomes at last day of follow-up - write to file 
	
	file write tablech ("Diabetes") _tab ("`a'") _tab ("Photoocagulation") _tab 
	* cumulative outcome - photocoagulation 
	sum _at1 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 
	* cumulative outcome - Antivegf 
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
	sum _at2 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 

	*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
	********************* Need to sort out graph axes
	twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
					(rarea _at2_lci _at2_uci days, color(blue%25)) ///
					(line _at1 days, sort lcolor(red)) ///
					(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
					, legend(order(1 "Photocoagulation" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
					title("Time to `a'", justification(left) size(med) )  	   ///
					yscale(range(0, 1)) 											///
					ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
					xlabel(0 (500) 2700, labsize(small))				   				///			
					ytitle("Cumulative outcomes (%)", size(medsmall)) ///
					xtitle("days since index date", size(medsmall))      		///
					graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

	graph export "$projdir/output/adjcurv_dm_`a'.svg", as(svg) replace

	* Close window 
	graph close
	drop days-date
}


use "$savedir\an_nodm_main_analysis", clear
logistic antivegf age_at_index i.gender i.eth5 i.imd i.smokstatus i.bmi_cat i.bl_af i.bl_depression i.bl_hypertension i.bl_kidney_failure i.bl_mi i.bl_neph_syndrome i.bl_stroke i.bl_hf i.bl_pad i.bl_copd i.bl_cancer yr_index i.index_statin i.index_acei i.index_antiplatelet i.index_arb i.index_betablocker i.index_ccb i.index_loop_diuretic i.index_mra i.index_oac i.index_otherantihypertensive i.index_ppi i.index_nsaid scr_meas_yr_prior tot_appts_yr_prior


* Estimate propensity scores
predict propensity

graph tw kdensity propensity if antivegf == 0 || kdensity propensity if antivegf == 1
*graph export "$projdir\output\ps_dm_all.jpg", replace
* Summarise distributions - 
bys antivegf: sum propensity, d
* Balance diagnostics
gen wts = 1/propensity if antivegf == 1
replace wts = 1/(1-propensity) if antivegf == 0

* Trim propensity score
egen min_ps_1 = min(propensity) if antivegf==1
egen min_ps_0 = min(propensity) if antivegf==0
egen max_ps_1 = max(propensity) if antivegf==1
egen max_ps_0 = max(propensity) if antivegf==0
foreach var in min_ps_1 min_ps_0 max_ps_1 max_ps_0 {
	sort `var'
	replace `var' = `var'[_n-1] if `var'==.
}

gen trim_min = max(min_ps_1, min_ps_0)
gen trim_max = min(max_ps_1, max_ps_0)

count if propensity<trim_min 
count if propensity>trim_max
* drop rather than generate new variable 
drop if (propensity<trim_min | propensity>trim_max)

* Calculate ATT weights 
gen att_weight = antivegf + (1-antivegf)*(propensity/(1-propensity))

* Without weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs
* With weighting
covbal antivegf age_at_index gender eth5 imd smokstatus bmi_cat bl_af bl_depression bl_hypertension bl_kidney_failure bl_mi bl_neph_syndrome bl_stroke bl_hf bl_pad bl_copd bl_cancer yr_index index_statin index_acei index_antiplatelet index_arb index_betablocker index_ccb index_loop_diuretic index_mra index_oac index_otherantihypertensive index_ppi index_nsaid scr_meas_yr_prior tot_appts_yr_prior, abs wt(att_weight)


* Cumulative hazard plots 
local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure event_mi event_neph_syndrome event_stroke event_hf event_pad event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure end_mi end_neph_syndrome end_stroke end_hf end_pad end_cvd_death end_zoster"
local n: word count `event'
forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'
	stset `b' [pweight=att_weight], failure(`a') origin(index_date) enter(index_date) id(patid)
	
	* Add cumulative incidence plots 
	* Setting df (degrees of freedom for restricted cubic splines) as 3 as this is default 
	* Setting dftvc (degrees of freedom for time-dependent effects) as 1 = linear effect of log time 
	stpm2 antivegf, dftvc(1) df(2) scale(hazard) eform
	summ _t
	local tmax=r(max)
	local tmaxplus1=r(max)+1

	range days 0 `tmax' `tmaxplus1'
	stpm2_standsurv if antivegf == 1, at1(antivegf 0) at2(antivegf 1) timevar(days) ci contrast(difference) fail

	gen date = index_date + days
	format date %tddd_Month

	for var _at1 _at2 _at1_lci _at1_uci _at2_lci _at2_uci _contrast2_1 _contrast2_1_lci _contrast2_1_uci: replace X=100*X
	
	*cumulative outcomes at last day of follow-up - write to file 
	
	file write tablech ("No diabetes") _tab ("`a'") _tab ("Cataract") _tab 
	* cumulative outcome - photocoagulation 
	sum _at1 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at1_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 
	* cumulative outcome - Antivegf 
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab  
	sum _at2 if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_lci if days==`tmax'
	file write tablech (r(mean)) _tab 
	sum _at2_uci if days==`tmax'
	file write tablech (r(mean)) _tab _n 

	*l date days_ph _at1 _at1_lci _at1_uci _at2 _at2_lci _at2_uci if days_ph<.
	********************* Need to sort out graph axes
	twoway  (rarea _at1_lci _at1_uci days, color(red%25)) ///
					(rarea _at2_lci _at2_uci days, color(blue%25)) ///
					(line _at1 days, sort lcolor(red)) ///
					(line _at2 days, sort lcolor(blue) lpattern(dash)) ///
					, legend(order(1 "Cataract" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
					title("Time to `a'", justification(left) size(med) )  	   ///
					yscale(range(0, 1)) 											///
					ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small))	///
					xlabel(0 (500) 2700, labsize(small))				   				///			
					ytitle("Cumulative outcomes (%)", size(medsmall)) ///
					xtitle("days since index date", size(medsmall))      		///
					graphregion(fcolor(white)) saving(adjcurv_`outcome', replace)

	graph export "$projdir/output/adjcurv_nodm_`a'.svg", as(svg) replace

	* Close window 
	graph close
	drop days-date
}

file close tablech