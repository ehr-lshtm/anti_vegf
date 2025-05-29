/*==============================================================================
DO FILE NAME:			3_forest_plots.do
DATE: 					20/02/2025
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	creates forest plots

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

import delimited "$projdir/output/primary_cox_models.txt", varnames(1) clear

keep if exposuregroup == "antivegf"

* log estimates 
gen hr_log = log(hr)
gen lci_log = log(lci)
gen uci_log = log(uci)

replace outcome = "40% reduction eGFR" if outcome=="egfr_40"
replace outcome = "Sustained 40% reduction eGFR" if outcome=="egfr_40_sustained"
replace outcome = "Progression of albuminuria stage" if outcome=="acr_increased"
replace outcome = "Atrial fibrillation" if outcome=="event_af"
replace outcome = "Kidney failure" if outcome=="event_kidney_failure_15"
replace outcome = "Myocardial infarction" if outcome=="event_mi"
replace outcome = "Nephrotic syndrome" if outcome=="event_neph_syndrome"
replace outcome = "Stroke" if outcome=="event_stroke"
replace outcome = "Heart failure" if outcome=="event_hf"
replace outcome = "Peripheral arterial disease" if outcome=="event_pad"
replace outcome = "Cardiovascular event" if outcome=="event_cv"
replace outcome = "Cardiovascular death" if outcome=="event_cvd_death"
replace outcome = "Zoster infection (negative control)" if outcome=="event_zoster"
replace outcome = "Incident hypertension" if outcome=="event_hypertension"
drop if outcome=="Stroke" | outcome=="Myocardial infarction" | outcome=="Peripheral arterial disease"
* order results 
* Re-order
gen order = 1 if outcome == "40% reduction eGFR" 
replace order = 2 if outcome == "Sustained 40% reduction eGFR" 
replace order = 3 if outcome == "Progression of albuminuria stage"
replace order = 6 if outcome =="Atrial fibrillation"
replace order = 5 if outcome =="Kidney failure" 
replace order = 9 if outcome =="Cardiovascular event" 
replace order = 4 if outcome =="Nephrotic syndrome" 
replace order = 8 if outcome =="Heart failure" 
replace order = 10 if outcome =="Cardiovascular death" 
replace order = 11 if outcome =="Zoster infection (negative control)" 
replace order = 7 if outcome =="Incident hypertension" 
sort group order

metan hr_log lci_log uci_log, eform  ///
	effect(Hazard Ratio) notable forestplot(null(1) dp(2) xlab(.5 1 2 3, force) favours("Favours Anti-VEGF                 "   #   "                 Favours control", nosymmetric) xtitle(, size(tiny)) graphregion(margin(zero) color(white)) texts(100) astext(65)) by(group) nowt nosubgroup nooverall nobox scheme(sj) label(namevar=outcome)  

graph export "$projdir/output/fp_primary_cox_models.svg", replace

* Stratified models 
import delimited "$projdir/output/primary_cox_models_stratified.txt", varnames(1) clear

* log estimates 
gen hr_log = log(hr)
gen lci_log = log(lci)
gen uci_log = log(uci)

replace outcome = "40% reduction eGFR" if outcome=="egfr_40"
replace outcome = "Sustained 40% reduction eGFR" if outcome=="egfr_40_sustained"
replace outcome = "Progression of albuminuria stage" if outcome=="acr_increased"
replace outcome = "Atrial fibrillation" if outcome=="event_af"
replace outcome = "Kidney failure" if outcome=="event_kidney_failure_15"
replace outcome = "Myocardial infarction" if outcome=="event_mi"
replace outcome = "Nephrotic syndrome" if outcome=="event_neph_syndrome"
replace outcome = "Stroke" if outcome=="event_stroke"
replace outcome = "Heart failure" if outcome=="event_hf"
replace outcome = "Peripheral arterial disease" if outcome=="event_pad"
replace outcome = "Cardiovascular event" if outcome=="event_cv"
replace outcome = "Cardiovascular death" if outcome=="event_cvd_death"
replace outcome = "Zoster infection (negative control)" if outcome=="event_zoster"
replace outcome = "Incident hypertension" if outcome=="event_hypertension"
drop if outcome=="Stroke" | outcome=="Myocardial infarction" | outcome=="Peripheral arterial disease"
 

replace strata = "eGFR stage 3a" if strata=="egfr mild"
replace strata = "eGFR stage 3b, 4 or 5" if strata=="egfr moderate"
replace strata = "Male" if strata == "gender 1"
replace strata = "Female" if strata == "gender 2"
replace strata = "Ethnicity: White" if strata == "eth_bin 0"
replace strata = "Ethnicity: Non-White" if strata == "eth_bin 1"
replace strata = "Ethnicity: Missing" if strata == "eth_bin 5"

gen strata_n = 0 
replace strata_n = 1 if strpos(strata, "eGFR")
replace strata_n = 2 if strpos(strata, "ACR")
replace strata_n = 3 if strpos(strata, "IMD")
replace strata_n = 4 if strpos(strata, "Eth")
replace strata_n = 5 if strpos(strata, "ale")

* Redact results with events <=5
gen redact_1 = events<=5 
bys group outcome strata_n : egen redact = max(redact_1)
replace hr_log = 0 if redact==1

keep if exposuregroup == "antivegf"

* order results 
* Re-order
gen order = 1 if outcome == "40% reduction eGFR" 
replace order = 2 if outcome == "Sustained 40% reduction eGFR" 
replace order = 3 if outcome == "Progression of albuminuria stage"
replace order = 6 if outcome =="Atrial fibrillation"
replace order = 5 if outcome =="Kidney failure" 
replace order = 9 if outcome =="Cardiovascular event" 
replace order = 4 if outcome =="Nephrotic syndrome" 
replace order = 8 if outcome =="Heart failure" 
replace order = 10 if outcome =="Cardiovascular death" 
replace order = 11 if outcome =="Zoster infection (negative control)" 
replace order = 7 if outcome =="Incident hypertension" 
sort group order strata

forvalues i=1/5 {
	metan hr_log lci_log uci_log if strata_n==`i' & group=="Diabetes", eform  ///
	effect(Hazard Ratio) notable forestplot(null(1) dp(2) xlab(.25 .5 1 2 3, force) favours("Favours Anti-VEGF             "   #   "             Favours control", nosymmetric) xtitle(, size(tiny)) graphregion(margin(zero) color(white)) texts(100) astext(65)) by(outcome) nowt nosubgroup nooverall nobox scheme(sj) label(namevar=strata)  
	graph export "$projdir/output/fp_stratified_cox_models_dm_`i'.svg", replace
}

forvalues i=1/5 {
	metan hr_log lci_log uci_log if strata_n==`i' & group=="No diabetes" & redact==0, eform  ///
	effect(Hazard Ratio) notable forestplot(null(1) dp(2) xlab(.25 .5 1 2 3, force) favours("Favours Anti-VEGF             "   #   "             Favours control", nosymmetric) xtitle(, size(tiny)) graphregion(margin(zero) color(white)) texts(100) astext(65)) by(outcome) nowt nosubgroup nooverall nobox scheme(sj) label(namevar=strata)  
	graph export "$projdir/output/fp_stratified_cox_models_nodm_`i'.svg", replace
}

* Sensitivity analyses 
foreach value in 1 3 4 6 {
	import delimited "$projdir/output/primary_cox_models_sensitivity_`value'.txt", varnames(1) clear

	keep if exposuregroup == "antivegf"

	* log estimates 
	gen hr_log = log(hr)
	gen lci_log = log(lci)
	gen uci_log = log(uci)

	replace outcome = "40% reduction eGFR" if outcome=="egfr_40"
replace outcome = "Sustained 40% reduction eGFR" if outcome=="egfr_40_sustained"
replace outcome = "Progression of albuminuria stage" if outcome=="acr_increased"
replace outcome = "Atrial fibrillation" if outcome=="event_af"
replace outcome = "Kidney failure" if outcome=="event_kidney_failure_15"
replace outcome = "Myocardial infarction" if outcome=="event_mi"
replace outcome = "Nephrotic syndrome" if outcome=="event_neph_syndrome"
replace outcome = "Stroke" if outcome=="event_stroke"
replace outcome = "Heart failure" if outcome=="event_hf"
replace outcome = "Peripheral arterial disease" if outcome=="event_pad"
replace outcome = "Cardiovascular event" if outcome=="event_cv"
replace outcome = "Cardiovascular death" if outcome=="event_cvd_death"
replace outcome = "Zoster infection (negative control)" if outcome=="event_zoster"
replace outcome = "Incident hypertension" if outcome=="event_hypertension"
drop if outcome=="Stroke" | outcome=="Myocardial infarction" | outcome=="Peripheral arterial disease"
* order results 
* Re-order
gen order = 1 if outcome == "40% reduction eGFR" 
replace order = 2 if outcome == "Sustained 40% reduction eGFR" 
replace order = 3 if outcome == "Progression of albuminuria stage"
replace order = 6 if outcome =="Atrial fibrillation"
replace order = 5 if outcome =="Kidney failure" 
replace order = 9 if outcome =="Cardiovascular event" 
replace order = 4 if outcome =="Nephrotic syndrome" 
replace order = 8 if outcome =="Heart failure" 
replace order = 10 if outcome =="Cardiovascular death" 
replace order = 11 if outcome =="Zoster infection (negative control)" 
replace order = 7 if outcome =="Incident hypertension" 
	sort group order

	* Redact results with events <=5
gen redact_1 = events<=5 
bys group outcome : egen redact = max(redact_1)
replace hr_log = 0 if redact==1
	
	metan hr_log lci_log uci_log, eform  ///
		effect(Hazard Ratio) notable forestplot(null(1) dp(2) xlab(.5 1 2 3, force) favours("Favours Anti-VEGF             "   #   "             Favours control", nosymmetric) xtitle(, size(tiny)) graphregion(margin(zero) color(white)) texts(100) astext(65)) by(group) nowt nosubgroup nooverall nobox scheme(sj) label(namevar=outcome)  

	graph export "$projdir/output/fp_primary_cox_models_sens_`value'.svg", replace
}

* Sensitivity analysis 5
import delimited "$projdir/output/primary_cox_models_sensitivity_5.txt", varnames(1) clear
keep if exposuregroup == "antivegf"
destring lci, replace

* log estimates 
gen hr_log = log(hr)
gen lci_log = log(lci)
gen uci_log = log(uci)

replace outcome = "40% reduction eGFR" if outcome=="egfr_40"
replace outcome = "Sustained 40% reduction eGFR" if outcome=="egfr_40_sustained"
replace outcome = "Progression of albuminuria stage" if outcome=="acr_increased"
replace outcome = "Atrial fibrillation" if outcome=="event_af"
replace outcome = "Kidney failure" if outcome=="event_kidney_failure_15"
replace outcome = "Myocardial infarction" if outcome=="event_mi"
replace outcome = "Nephrotic syndrome" if outcome=="event_neph_syndrome"
replace outcome = "Stroke" if outcome=="event_stroke"
replace outcome = "Heart failure" if outcome=="event_hf"
replace outcome = "Peripheral arterial disease" if outcome=="event_pad"
replace outcome = "Cardiovascular event" if outcome=="event_cv"
replace outcome = "Cardiovascular death" if outcome=="event_cvd_death"
replace outcome = "Zoster infection (negative control)" if outcome=="event_zoster"
replace outcome = "Incident hypertension" if outcome=="event_hypertension"
drop if outcome=="Stroke" | outcome=="Myocardial infarction" | outcome=="Peripheral arterial disease"
* order results 
* Re-order
gen order = 1 if outcome == "40% reduction eGFR" 
replace order = 2 if outcome == "Sustained 40% reduction eGFR" 
replace order = 3 if outcome == "Progression of albuminuria stage"
replace order = 6 if outcome =="Atrial fibrillation"
replace order = 5 if outcome =="Kidney failure" 
replace order = 9 if outcome =="Cardiovascular event" 
replace order = 4 if outcome =="Nephrotic syndrome" 
replace order = 8 if outcome =="Heart failure" 
replace order = 10 if outcome =="Cardiovascular death" 
replace order = 11 if outcome =="Zoster infection (negative control)" 
replace order = 7 if outcome =="Incident hypertension" 
sort group order

gen time2 = floor(time)

foreach time in 183 365 730 1095 1461 1826 {
	metan hr_log lci_log uci_log if time2==`time', eform  effect(Hazard Ratio) notable forestplot(null(1) dp(2) xlab(.25 .5 1 2 3, force) favours("Favours Anti-VEGF             "   #   "             Favours control", nosymmetric) xtitle(, size(tiny)) graphregion(margin(zero) color(white)) texts(100) astext(65)) by(group) nowt nosubgroup nooverall nobox scheme(sj) label(namevar=outcome)  
	graph export "$projdir/output/fp_cox_models_sens_5_`time'.svg", replace
}
