/*==============================================================================
DO FILE NAME:			2_define_outcomes.do
DATE: 					03/12/2024
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	defines outcomes and generates analysis files 

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\2_define_outcomes_$S_DATE.log", append

* Generate closest egfr to index date
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\eGFR", clear
	duplicates drop patid obsdate SCr unit, force
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) keep(match) nogen
	merge m:1 patid using "$savedir\\`grp'\cr_bl_egfr"
	* Check if there are multiple measures on the same day
	bys patid obsdate: gen total = _N
	tab total 
	bys patid obsdate: egen min_egfr = min(egfr)
	bys patid obsdate: egen max_egfr = max(egfr)
	gen diff_egfr = min_egfr - max_egfr	
	gen outside_fu = (obsdate <= index_date | obsdate>end_fu)
	drop if outside_fu
	drop outside_fu
		
	* Rules for mulitple measures on same day - none currently
	gen egfr_percent_change = ((egfr-prior_egfr)/prior_egfr)*100
	* Flag >=40% reduction
	gen egfr_40 = egfr_percent_change <= -40 
	* Identify sustained reduction
	bys patid (obsdate): gen total_changes = sum(egfr_40)
	sort patid obsdate 
	gen egfr_40_sustained = total_changes>total_changes[_n-1] & patid==patid[_n-1] & total_changes[_n-1]>total_changes[_n-2] & patid==patid[_n-2]
	save "$savedir\\`grp'\all_egfr_fu", replace
	* Need to generate time-varying files with end date of first 40% reduction then another file with first sustained 40% reduction 
	* 40% reduction end date
	keep if egfr_40==1
	* 40% reduction end date
	bys patid: egen end_egfr_40 = min(obsdate)
	bys patid: egen end_egfr_40_sustained = min(obsdate) if egfr_40_sustained==1
	
	keep patid egfr_40* end* 
	duplicates drop 
	* Doesn't take out all duplicates as rows where sustained equals 1 and 0
	bys patid: egen any_sustained = max(egfr_40_sustained)
	drop if egfr_40_sustained==0 & any_sustained==1
	count
	codebook patid 
	merge 1:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) nogen
	replace end_egfr_40 = end_fu if end_egfr_40==.
	replace end_egfr_40_sustained = end_fu if end_egfr_40_sustained==.
	replace egfr_40 = 0 if egfr_40==.
	replace egfr_40_sustained = 0 if egfr_40_sustained==. 
	save "$savedir\\`grp'\outcome_egfr_40", replace
}
	
* Progression of albuminuria stage
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\ACR", clear
	duplicates drop patid obsdate ACR unit, force
	merge m:1 patid using "$savedir\\`grp'\cr_bl_acr", nogen
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) keep(match) nogen
	* Define albuminuria categories
	egen acr_cat= cut(ACR), at(0, 30, 300, 5000) icodes
	*label define ACR 0"A1 Normal to mildly increased" 1 "A2 Moderately increased" 2 "A3 Severely increased"
	label values prior_acr_cat ACR
	label var acr_cat "Albuminuria category"
	* Keep only those measures within follow-up
	gen outside_fu = (obsdate <= index_date | obsdate>end_fu)
	drop if outside_fu
	drop outside_fu	
	
	* Check if there are multiple measures on the same day
	bys patid obsdate: gen total = _N
	tab total 
	bys patid obsdate: egen min_acr = min(acr_cat)
	bys patid obsdate: egen max_acr = max(acr_cat)
	gen diff_acr = min_acr - max_acr
	tab diff_acr
	* Only a few where stage is different between measures on the same day, will take lowest stage on day 
	gen drop = (acr_cat!=min_acr & total!=1)
	drop if drop 
	duplicates drop patid obsdate acr_cat, force 
	drop total 
	bys patid obsdate: gen total = _N
	tab total
	* Determine if changed category since the last measure
	sort patid obsdate 
	bys patid (obsdate): gen acr_increased = acr_cat[_n-1] < acr_cat
	bys patid (obsdate): replace acr_increased=1 if prior_acr_cat < acr_cat & _n==1
	* Set to missing if no prior measure 
	bys patid (obsdate): replace acr_increased=. if prior_acr_cat==. & _n==1
	* Set end dates 
	keep if acr_increased==1
	bys patid: egen end_acr = min(obsdate)
	
	keep patid acr_increased end* 
	duplicates drop 
	count
	codebook patid 
	merge 1:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) nogen
	replace end_acr = end_fu if end_acr==.
	replace acr_increased = 0 if acr_increased==.
	save "$savedir\\`grp'\outcome_acr", replace

}
	
* 	New kidney failure, hypertension and CV events 
foreach grp in antivegf cataract photocoag {
	
	display in red "*******************Observation file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(patid yob index_date end_fu) keep(match) nogen
	
	foreach script in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad {
	
		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_`script'", nogen keepusing(`script'_cprd)
	
	}
	egen keep = rowmax(af_cprd-pad_cprd)
	tab keep 
	keep if keep==1
	
	save `tempfile', replace
	
		/*******************************************************************************
	#A2. Loop through subsequent (from 2 onwards) separate test extract files in 
		turn and append the results to the first extract file saved in #1
	*******************************************************************************/
	if "`grp'"=="antivegf" {
		local j 5
	}
	if "`grp'"=="cataract" {
		local j 74
	}
	if "`grp'"=="photocoag" {
		local j 3
	}
	di `j'
	forvalues n=2/`j' {
		display in red "*******************Observation file number: `n'*******************"

		use "$rawdata\\`grp'\observation_`n'", clear
		
		merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(patid yob index_date end_fu) keep(match) nogen
			
		foreach script in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad  {
	
		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_`script'", nogen keepusing(`script'_cprd)
	
		}
		egen keep = rowmax(af_cprd-pad_cprd)
		tab keep 
		keep if keep==1

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}

	
	gen eventdate = date(obsdate, "DMY")
	* Drop missing event dates
	drop if eventdate==.
	* Keep records after index and before end of follow-up
	gen after_index = eventdate>index_date & eventdate<=end_fu
	tab after_index 
	keep if after_index
	
	* Spread vars so single row per person
	foreach script in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad  {
			bys patid: egen event_`script' = max(`script'_cprd)
			bys patid: egen first_code_`script' = min(eventdate) if `script'_cprd==1
			format first_code_`script' %dD/N/CY
			* Spread to all rows
			bys patid (first_code_`script'): replace first_code_`script'=first_code_`script'[1] if first_code_`script'==.
		}
	
	keep patid event_af-first_code_pad
	duplicates drop 
	codebook patid 
	merge 1:1 patid using "$savedir\\`grp'\cr_vars", keepusing(patid) nogen
	* Tidying 
	foreach script in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad  {
		replace event_`script'=0 if event_`script'==.
	}
	codebook patid 
	save "$savedir\\`grp'\cr_aurum_events", replace
	
}
	
* Find hospital events
use "$savedir/all_hosp_bl", clear
foreach grp in antivegf cataract photocoag {
	preserve 
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) keep(match) nogen
	keep if (admitdate>index_date & admitdate<=end_fu)
	foreach var in mi stroke hf arrhythmia kidney_failure zoster {
		* Determine if row is prior to index 
		*gen event_hosp_`var'_i = (admitdate>index_date & admitdate<=end_fu & hes_`var'==1)
		* Drop records after index 
		*drop if (event_hosp_`var'_i==0 & hes_`var'!=0)
		* Create static flag for prior comorbidities
		bys patid: egen event_hosp_`var' = max(hes_`var')
		* Identify first code for each comorbidity
		bys patid: egen first_code_`var' = min(admitdate) if hes_`var'==1
		bys patid: egen first_code_`var'_hes = max(first_code_`var')
		drop first_code_`var'
	}
	count if admitdate<=index_date
	drop if event_hosp_mi==0 & event_hosp_stroke==0 & event_hosp_hf==0 & event_hosp_arrhythmia==0 & event_hosp_kidney_failure==0 & event_hosp_zoster==0
	keep patid event_hosp_mi event_hosp_stroke event_hosp_hf event_hosp_arrhythmia event_hosp_kidney_failure event_hosp_zoster first_code*
	duplicates drop 
	count 
	codebook patid 
	save "$savedir\\`grp'\cr_hes_events", replace
	restore 
}

* Merge together with comorbidities in GP records
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\cr_aurum_events", clear
	merge 1:1 patid using "$savedir\\`grp'\cr_hes_events"
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(end_fu) keep(match) nogen
	* Edit some variable names to match Aurum variables 
	rename event_hosp_arrhythmia event_hosp_af 
	rename first_code_arrhythmia_hes first_code_af_hes
	foreach var in mi stroke hf af kidney_failure zoster {
		replace event_`var' = event_hosp_`var' if event_`var'==0 & event_hosp_`var'==1
		replace first_code_`var' = first_code_`var'_hes if first_code_`var'_hes<first_code_`var'
		drop event_hosp_`var' first_code_`var'_hes 
		count if event_`var'==1 & first_code_`var'==.
	}
	foreach var in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad {
		rename first_code_`var' end_`var'
	}
	save "$savedir\\`grp'\cr_events", replace
}	

* Determine if egfr <15 during follow-up 
* Generate closest egfr to index date
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\eGFR", clear
	duplicates drop patid obsdate SCr unit, force
	merge m:1 patid using "$savedir\\`grp'\cr_vars", keepusing(index_date end_fu) keep(match) nogen
	* Check if there are multiple measures on the same day
	bys patid obsdate: gen total = _N
	tab total 
	bys patid obsdate: egen min_egfr = min(egfr)
	bys patid obsdate: egen max_egfr = max(egfr)
	gen diff_egfr = min_egfr - max_egfr	
	gen outside_fu = (obsdate <= index_date | obsdate>end_fu)
	drop if outside_fu
	drop outside_fu
	gen egfr_15 = egfr<15 
	tab egfr_15 
	keep if egfr_15 
	bys patid: egen first_egfr_15 = min(obsdate)
	keep if first_egfr_15 == obsdate
	keep patid obsdate egfr_15 first_egfr_15 
	duplicates drop 
	count 
	codebook patid 
	drop obsdate
	save "$savedir\\`grp'\cr_egfr_15", replace
}
	
* CV death 
use  "$rawdata\linked_data_20250106\linked_death", clear 
merge m:1 patid using "$savedir\cr_all_study_pop", keep(match) keepusing(regenddate cprd_ddate) nogen
* File has multiple rows per patient. Only need underlying cause of death and date of death for now will deduplicate based on death date - choosing earliest date if more than one date
drop if s_underlying_cod_icd10==""
bys patid: gen n = _N
tab n 
* Format dates 
gen date_of_death = date(reg_date_of_death, "YMD")
gen cprd_death_date = date(cprd_ddate, "DMY")
format date_of_death cprd_death_date %dD/N/CY
bys patid (date_of_death): gen diff_ons = date_of_death!=date_of_death[_n-1] & patid==patid[_n-1]
* Determine first death date
bys patid: egen death_date = min(date_of_death)
* Only two people who haev multiple ONS death dates and one of those dates matches CPRD death date - in both cases only one day out from minimum therefore keeping as minimum  
* Identify CV underlying causes
gen icd_letter = substr(s_underlying_cod_icd10, 1, 1)
gen icd_number = substr(s_underlying_cod_icd10, 2, 3)
gen icd_number_2 = icd_number + substr("000", 1, 3 - length(icd_number))
gen event_cvd_death = 1 if icd_letter=="I"
keep if event_cvd_death==1
rename death_date end_cvd_death
keep patid event_cvd_death end_cvd_death 
duplicates drop 
codebook patid 
count 
drop if end_cvd_death > date("01Mar2020", "DMY")
save "$savedir\cr_cvd_death_dates", replace 
	
* Add outcomes to baseline dataset
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\cr_vars", clear
	merge 1:1 patid using "$savedir\\`grp'\cr_events", nogen 
	drop _merge
	merge 1:1 patid using "$savedir\cr_cvd_death_dates"
	drop if _merge==2
	drop _merge 
	replace event_cvd_death = 0 if end_cvd_death > end_fu & end_cvd_death!=.
	replace end_cvd_death = . if end_cvd_death > end_fu & end_cvd_death!=.
	foreach script in af hypertension kidney_failure mi neph_syndrome stroke zoster hf pad cvd_death {
		replace event_`script'=0 if event_`script'==.
		replace end_`script' = end_fu if end_`script'==.
	}
	merge 1:1 patid using "$savedir\\`grp'\outcome_acr", nogen 
	merge 1:1 patid using "$savedir\\`grp'\outcome_egfr_40", nogen keep(match)
	merge 1:1 patid using "$savedir\\`grp'\cr_egfr_15", nogen 
	gen event_kidney_failure_15 = event_kidney_failure 
	replace event_kidney_failure_15 = 1 if egfr_15==1 & first_egfr_15 < end_kidney_failure
	gen end_kidney_failure_15 = end_kidney_failure 
	replace end_kidney_failure_15 = first_egfr_15 if egfr_15==1 & first_egfr_15 < end_kidney_failure
	save "$savedir\\`grp'\cr_bl_plus_outcomes", replace
}

* Create cohorts for analysis
use "$savedir\antivegf\cr_bl_plus_outcomes", clear 
keep if bl_dm==1 
gen antivegf=1
append using "$savedir\photocoag\cr_bl_plus_outcomes"
replace antivegf = 0 if antivegf==.
* Update cohort based on linked HES data received after defining cohort as CPRD did not identify procedures that occurred outside of study period
merge 1:1 patid using "$savedir\antivegf\hes_first_code", gen(merge_antivegf)
drop if merge_antivegf==2 
merge 1:1 patid using "$savedir\photocoag\hes_first_code", gen(merge_photocoag)
drop if merge_photocoag==2 
gen out = ((antivegf==1 & merge_photocoag==3 & first_photocoag < index_date) | (antivegf==0 & merge_antivegf==3 & first_antivegf < index_date) | (antivegf==1 & first_antivegf==.) | (antivegf==0 & first_photocoag==.))
tab antivegf out
* Generate composite CV event outcome 
egen event_cv = rowmax(event_mi event_stroke event_pad)
egen end_cv = rowmin(end_mi end_stroke end_pad)
save "$savedir\an_dm_main_analysis", replace

use "$savedir\antivegf\cr_bl_plus_outcomes", clear 
keep if bl_dm==0 
gen antivegf=1
append using "$savedir\cataract\cr_bl_plus_outcomes"
replace antivegf = 0 if antivegf==.
* Update cohort based on linked HES data received after defining cohort as CPRD did not identify procedures that occurred outside of study period
merge 1:1 patid using "$savedir\antivegf\hes_first_code", gen(merge_antivegf)
drop if merge_antivegf==2 
merge 1:1 patid using "$savedir\cataract\hes_first_code", gen(merge_cataract)
drop if merge_cataract==2 
gen out = ((antivegf==1 & merge_cataract==3 & first_cataract < index_date) | (antivegf==0 & merge_antivegf==3 & first_antivegf < index_date) | (antivegf==1 & first_antivegf==.) | (antivegf==0 & first_cataract==.))
tab antivegf out
* Generate composite CV event outcome 
egen event_cv = rowmax(event_mi event_stroke event_pad)
egen end_cv = rowmin(end_mi end_stroke end_pad)
save "$savedir\an_nodm_main_analysis", replace