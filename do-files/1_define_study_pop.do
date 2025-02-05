/*==============================================================================
DO FILE NAME:			1_define_study_pop.do
DATE: 					16/10/2024
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	defines study populations

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

* run globals.do first

log using "$projdir\logs\1_refine_$S_DATE.log", append

/*Source populations criteria applied to all 3 cohorts:
	1. age 18+, 2. acceptable, 3. >=1 year registration, 4. >=1 eGFR
	*/

* Will need to know whether people have had prior eye disease treatment so merge together three files with first code 
use "$rawdata\antivegf\vegf_first_date", clear 
rename first_code first_code_antivegf 
merge 1:1 patid using "$rawdata\photocoag\photocoag_first_date", nogen
rename first_code first_code_photocoag
merge 1:1 patid using "$rawdata\cataract\cataract_first_date", nogen
rename first_code first_code_cataract
* Identify first treatment across 3 types 
gen first_rx = first_code_antivegf
replace first_rx = first_code_photocoag if first_code_photocoag<first_code_antivegf
replace first_rx = first_code_cataract if first_code_cataract<first_code_photocoag & first_code_cataract<first_code_antivegf
* Flags for inclusion 
foreach grp in antivegf photocoag cataract {
	gen inclusion_prior_rx_`grp' = (first_code_`grp' == first_rx)
}
* Small number meet inclusion for multiple cohorts (more than one treatment on same day) set these not not be included
egen multi = rowtotal( inclusion_prior_rx_antivegf- inclusion_prior_rx_cataract)
tab multi 
foreach grp in antivegf photocoag cataract {
	replace inclusion_prior_rx_`grp' = 0 if multi>1
}
save "$savedir\first_code_inclusions", replace 
	
* First merge patient and practice file variables to determine 1st 3 criteria
mkdir "$savedir/antivegf"

use "$rawdata\antivegf\patient_1", clear 

* merge on first injection dates
merge 1:1 patid using "$rawdata\antivegf\vegf_first_date", gen(match_code)
tab match_code 
keep if match_code==3
drop match_code 
* Take out those who were no longer registered with their practice at first code
gen reg_end = date(regenddate, "DMY")
gen out = (reg_end<=first_code)
tab out
drop if out==1
drop out
*1. age 18+
gen yr_index = year(first_code)
gen age = yr_index - yob 
sum age
gen inclusion_age = age>=18
* Take out people with indeterminate sex
replace inclusion_age=0 if gender==3
replace inclusion_age = 0 if inclusion_age==.

* 2. acceptable 
gen inclusion_accept = (acceptable == 1)
count 
replace inclusion_accept = 0 if inclusion_accept==.

* 3. >=1 year registration at first injection 
* 
* Format regstartdate 
gen reg_start = date(regstartdate, "DMY")
gen reg_yr = year(reg_start)
gen reg_mth = month(reg_start)
gen reg_day = day(reg_start)
gen reg_start_plus = dmy(reg_day, reg_mth, reg_yr+1)
* Need to update for leap years as these are missing 
replace reg_start_plus = dmy(reg_day-1, reg_mth, reg_yr+1) if reg_start_plus==.
format reg_start reg_start_plus %dD/N/CY

* Time since registration at first injection 
gen time_reg = first_code - reg_start_plus
sum time_reg
gen inclusion_reg = time_reg>=0
tab inclusion_age 
tab inclusion_accept
* save version of patients meeting inclusion criteria so far
save "$savedir\antivegf\prelim_pts", replace

**********************************************************************************
*********************COMPARATORS**************************************************

* Cataracts - no diabetes & Photocoagulation - diabetes 
*mkdir "$savedir/cataract" 
*mkdir "$savedir/photocoag" 
* Determine if have diabetes at first code
foreach grp in cataract photocoag {
	
	use "$rawdata\\`grp'\patient_1", clear 
	
	* merge on first injection dates and keep only matches
	merge 1:1 patid using "$rawdata\\`grp'\\`grp'_first_date", gen(match_code)
	tab match_code 
	keep if match_code==3
	drop match_code 
	gen reg_end = date(regenddate, "DMY")
	gen out = (reg_end<=first_code)
	tab out
	drop if out==1
	drop out
	
	save "$savedir\\`grp'\prelim_pts_all", replace

	**** Determine if have diabetes at first code
	tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	describe 
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_all", keepusing(patid gender yob first_code) keep(match) nogen

	merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_diabetes", gen(dm_code) keep(match)

	save `tempfile', replace

		/*******************************************************************************
	#A2. Loop through subsequent (from 2 onwards) separate observation extract files in 
		turn and append the results to the first extract file saved in #1
	*******************************************************************************/
	* Different number of files depending on which group therefore loop
	if "`grp'"=="cataract" {
		local j 74
	}
	else {
		local j 3
	}
	di `j'
	forvalues n=2/`j' {
		display in red "*******************Observation file number: `n'*******************"
		use "$rawdata\\`grp'\observation_`n'", clear
		
		merge m:1 patid using "$savedir\\`grp'\prelim_pts_all", keepusing(patid gender yob first_code) keep(match) nogen

		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_diabetes", gen(dm_code) keep(match)
			
		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}
	keep if bl_dm==1 
	save "$savedir\\`grp'\dm", replace
}

* Determine diabetes status and other inlcusion flags, for cataract cohort keep those without diabetes, for photocoag keep only those with diabetes 
foreach grp in cataract photocoag {
	use "$savedir\\`grp'\dm", clear
	keep patid obsdate bl_dm first_code
	duplicates drop
	tab bl_dm, m
	gen obsdate1=date(obsdate,"DMY")
	format obsdate1 %dD/N/CY

	drop obsdate 
	rename obsdate1 obsdate
	drop if obsdate==.
	* Identify diabetes prior to first injection 
	gen prior_first_code = obsdate<=first_code 
	tab prior_first_code 
	keep if prior_first_code
	drop obsdate prior_first_code first_code 
	duplicates drop 
	* Merge on DM drugs in 6 months prior to index
	merge 1:1 patid using "$savedir\\`grp'\prelim_pts_all"
	replace bl_dm=0 if bl_dm==.
	tab _merge bl_dm
	drop _merge
	
	tab bl_dm 
	if "`grp'"=="cataract" {
		drop if bl_dm==1
	}
	else if "`grp'"=="photocoag" {
		drop if bl_dm==0
	}
	
	
	
	*1. age 18+
	gen yr_index = year(first_code)
	gen age = yr_index - yob 
	sum age
	gen inclusion_age = age>=18
	* Take out people with indeterminate sex
	replace inclusion_age=0 if gender==3
	replace inclusion_age = 0 if inclusion_age==.

	* 2. acceptable 
	gen inclusion_accept = (acceptable == 1)
	count 
	replace inclusion_accept = 0 if inclusion_accept==.

	* 3. >=1 year registration at first injection 
	* 
	* Format regstartdate 
	gen reg_start = date(regstartdate, "DMY")
	gen reg_yr = year(reg_start)
	gen reg_mth = month(reg_start)
	gen reg_day = day(reg_start)
	gen reg_start_plus = dmy(reg_day, reg_mth, reg_yr+1)
	* Need to update for leap years as these are missing 
	replace reg_start_plus = dmy(reg_day-1, reg_mth, reg_yr+1) if reg_start_plus==.
	format reg_start reg_start_plus %dD/N/CY

	* Time since registration at first injection 
	gen time_reg = first_code - reg_start_plus
	sum time_reg
	gen inclusion_reg = time_reg>=0
	tab inclusion_age 
	tab inclusion_accept

* save version of patients meeting inclusion criteria so far
save "$savedir\\`grp'\prelim_pts", replace
}

* Get eGFR values for all three groups 
do "$dodir/1_prog_getSCr_Aurum.do"

* Identify eGFR values before first code  
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\eGFR", clear
	count 
	* Merge on egfr results
	merge m:1 patid using "$savedir\\`grp'\prelim_pts", keepusing(first_code) keep(match)
	count
	* Identify tests prior to first injection 
	gen prior_first_code = obsdate<=first_code 
	drop if prior_first_code==0 & _merge==3
	* Determine test closest to first injection
	bys patid: egen closest_egfr = max(obsdate)

	gen time_closest_egfr = first_code - closest_egfr 
	sum time_closest_egfr if prior_first_code==1
	format closest_egfr %dD/N/CY
	* Identify those who meet the inclusion criteria i.e. egfr within 18 months prior to first injection 
	gen inclusion_egfr = time_closest_egfr <=548
	replace inclusion_egfr = 0 if inclusion_egfr==.
	codebook patid 
	merge m:1 patid using "$savedir\\`grp'\prelim_pts", keepusing(patid yob gender inclusion*) gen(prelim) 
	codebook patid
	replace inclusion_egfr = 0 if prelim==2
	keep patid inclusion* yob gender first_code 
	duplicates drop
	count
	save "$savedir\\`grp'\prelim_pts_2", replace 
}

**** Determine if MI or stroke in 6 months prior to first code
foreach grp in antivegf cataract photocoag {
	/*tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	describe 
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_2", keepusing(patid gender yob first_code) keep(match) nogen

	merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_mi", gen(mi_code) keepusing(mi)

	merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_stroke", gen(stroke_code) keepusing(stroke)

	foreach var in mi_code stroke_code {
		drop if `var'==2
		*drop `var'
	}
	keep if mi_cprd==1 | stroke_cprd==1
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
		if "`grp'"=="photocoag"  {
			local j 3
		}
	di `j'
	forvalues n=2/`j' {
	display in red "*******************Observation file number: `n'*******************"
		use "$rawdata\\`grp'\observation_`n'", clear
		
		merge m:1 patid using "$savedir\\`grp'\prelim_pts_2", keepusing(patid first_code) keep(match) nogen

		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_mi", gen(mi_code) keepusing(mi)

		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_stroke", gen(stroke_code) keepusing(stroke)

		foreach var in mi_code stroke_code {
			drop if `var'==2
			*drop `var'
		}
		keep if mi_cprd==1 | stroke_cprd==1
		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}


	save "$savedir\\`grp'\mi_stroke", replace*/

	use "$savedir\\`grp'\mi_stroke", clear 
	drop yob gender
	duplicates drop

	* drop if eventdate missing 
	* but check if sysdate available and replace missing eventdate with sysdate if available
	replace obsdate=enterdate if (obsdate=="" & enterdate!="")
	gen obsdate1=date(obsdate,"DMY")
	format obsdate1 %dD/N/CY

	drop obsdate 
	rename obsdate1 obsdate
	drop if obsdate==.

	* Identify tests prior to first injection 
	gen prior_first_code = obsdate<=first_code 
	keep if prior_first_code

	* Determine test closest to first injection
	bys patid: egen closest_event = max(obsdate)

	gen time_closest_event = first_code - closest_event 
	sum time_closest_event if prior_first_code==1
	format closest_event %dD/N/CY
	* Identify those who meet the inclusion criteria i.e. MI/stroke with 6 months prior to first injection 
	gen inclusion_mi_stroke = time_closest_event >=182
	codebook patid if inclusion_mi_stroke==0
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_2", gen(pts) keepusing(patid yob gender first_code inclusion*)
	replace inclusion_mi_stroke=1 if pts==2
	codebook patid if inclusion_mi_stroke==1
	codebook patid if inclusion_mi_stroke==0
	keep patid inclusion* yob gender first_code
	duplicates drop
	merge 1:1 patid using "$savedir\first_code_inclusions", keep(match) keepusing(inclusion_prior_rx_`grp') nogen 
	save "$savedir\\`grp'\prelim_pts_3", replace 
}

**** Determine if have diabetes at first code for antivegf cohort 
tempfile tempfile 
use "$rawdata\antivegf\observation_1", clear
describe 
merge m:1 patid using "$savedir/antivegf/prelim_pts", keepusing(patid gender yob first_code) keep(match) nogen

merge m:1 medcodeid using "$savedir/codelists/cprd_aurum_codes_diabetes", gen(dm_code)

drop if dm_code==2

save `tempfile', replace

	/*******************************************************************************
#A2. Loop through subsequent (from 2 onwards) separate test extract files in 
	turn and append the results to the first extract file saved in #1
*******************************************************************************/

forvalues n=2/5 {
	display in red "*******************Observation file number: `n'*******************"
	use "$rawdata\antivegf\observation_`n'", clear
	
	merge m:1 patid using "$savedir/antivegf/prelim_pts", keepusing(patid gender yob first_code) keep(match) nogen

merge m:1 medcodeid using "$savedir/codelists/cprd_aurum_codes_diabetes", gen(dm_code)

drop if dm_code==2

	* add the file containing records for the specified comorbidity
	* to make one file containing all specified comorbidiy records for the
	* clinical extract specified
	append using `tempfile'
	
	* save
	save `tempfile', replace
}

keep if bl_dm==1 
save "$savedir\antivegf\dm", replace

duplicates drop

* drop if eventdate missing 
* but check if sysdate available and replace missing eventdate with sysdate if available
replace obsdate=enterdate if (obsdate=="" & enterdate!="")
gen obsdate1=date(obsdate,"DMY")
format obsdate1 %dD/N/CY

drop obsdate 
rename obsdate1 obsdate
drop if obsdate==.

* Identify tests prior to first injection 
gen prior_first_code = obsdate<=first_code 
keep if prior_first_code
keep patid bl_dm
duplicates drop 

codebook patid 
count 
merge m:1 patid using "$savedir\antivegf\prelim_pts_3", gen(pts) keepusing(patid yob gender first_code inclusion*)
replace bl_dm=0 if pts==2

keep patid inclusion* yob gender first_code bl_dm
duplicates drop
save "$savedir/antivegf/prelim_pts_3", replace 

* Determine if have MI or stroke in 6 months prior to index using linked data 
* All cohorts together in linked file
use "$rawdata\linked_data_20250106\linked_hes_diagnosis_hosp", clear 
merge m:1 icd using "$savedir/codelists/hes_mi_codes", gen(merge_mi_codes)
merge m:1 icd using "$savedir/codelists/hes_stroke_codes", gen(merge_stroke_codes)
keep if (merge_mi_codes==3 | merge_stroke_codes==3)
* Format date 
gen admitdate = date(admidate, "YMD")
format admitdate %dD/N/CY
drop admidate 
* Find events for each cohort
foreach grp in antivegf cataract photocoag {
	preserve 
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_3", keep(match) keepusing(first_code)
	* Identify tests prior to first injection 
	gen prior_first_code = admitdate<=first_code 
	keep if prior_first_code
	* Determine test closest to first injection
	bys patid: egen closest_event = max(admitdate)
	gen time_closest_event = first_code - closest_event 
	sum time_closest_event if prior_first_code==1
	format closest_event %dD/N/CY
	* Identify those who meet the inclusion criteria i.e. MI/stroke with 6 months prior to first injection 
	gen hes_mi_stroke = time_closest_event <182
	keep patid hes_mi_stroke 
	duplicates drop 
	codebook patid 
	count 
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_3", nogen
	replace inclusion_mi_stroke=0 if hes_mi_stroke==1
	drop hes_mi_stroke
	save "$savedir\\`grp'\prelim_pts_4", replace 
	restore 
}

* Check that none of patients died before first code using linked ONS data 
use  "$rawdata\linked_data_20250106\linked_death", clear 
* File has multiple rows per patient. Only need underlying cause of death and date of death for now will deduplicate based on death date - choosing earliest date if more than one date
drop if s_underlying_cod_icd10==""
bys patid: gen n = _N
tab n 
* Format date 
gen date_of_death = date(reg_date_of_death, "YMD")
format date_of_death %dD/N/CY
* Determine first death date
bys patid: egen death_date = min(date_of_death)
gen days = date_of_death-death_date
sum days if days!=0
keep patid death_date 
duplicates drop 
codebook patid 
count 
* merge each cohort file on to check did not die before first code 
foreach grp in antivegf cataract photocoag {
	preserve
	merge 1:1 patid using "$savedir\\`grp'\prelim_pts_4", keep(match) keepusing(first_code)
	
	gen inclusion_died_prior = death_date > first_code 
	tab inclusion_died_prior
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_4", nogen
	replace inclusion_died_prior=1 if inclusion_died_prior==.
	save "$savedir\\`grp'\prelim_pts_4", replace
	restore 
}

* Prelim flowchart numbers 
* Antivegf
use "$savedir/antivegf/prelim_pts_4", clear 
* DM
count if bl_dm==1
count if inclusion_died_prior==0 & bl_dm==1
count if inclusion_prior_rx_antivegf==0 & inclusion_died_prior==1 & bl_dm==1
count if inclusion_age==0 & inclusion_prior_rx_antivegf==1 & bl_dm==1
count if inclusion_accept==0 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==1
count if inclusion_reg==0 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==1
count if inclusion_egfr==0 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==1
count if inclusion_mi_stroke==0 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1  & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==1

count if inclusion_mi_stroke==1 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1  & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==1

* No DM
count if bl_dm==0
count if inclusion_died_prior==0 & bl_dm==0
count if inclusion_prior_rx_antivegf==0 & inclusion_died_prior==1 & bl_dm==0
count if inclusion_age==0 & inclusion_prior_rx_antivegf==1 & bl_dm==0
count if inclusion_accept==0 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==0
count if inclusion_reg==0 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==0
count if inclusion_egfr==0 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==0
count if inclusion_mi_stroke==0 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1  & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==0

count if inclusion_mi_stroke==1 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1  & inclusion_prior_rx_antivegf==1 & inclusion_died_prior==1 & bl_dm==0

* cataract 
use "$savedir/cataract/prelim_pts_4", clear 
count 
count if inclusion_died_prior==0
count if inclusion_prior_rx_cataract==0 & inclusion_died_prior==1
count if inclusion_age==0 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1
count if inclusion_accept==0 & inclusion_age==1 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1
count if inclusion_reg==0 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1
count if inclusion_egfr==0 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1
count if inclusion_mi_stroke==0 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1

count if inclusion_mi_stroke==1 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_cataract==1 & inclusion_died_prior==1

* Photocoagulation cohort 
use "$savedir/photocoag/prelim_pts_4", clear 
count 
count if inclusion_died_prior==0
count if inclusion_prior_rx_photocoag==0 & inclusion_died_prior==1
count if inclusion_age==0 & inclusion_prior_rx_photocoag==1 & inclusion_died_prior==1
count if inclusion_accept==0 & inclusion_age==1 & inclusion_prior_rx_photocoag==1 & inclusion_died_prior==1
count if inclusion_reg==0 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_photocoag==1  & inclusion_died_prior==1
count if inclusion_egfr==0 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_photocoag==1 & inclusion_died_prior==1
count if inclusion_mi_stroke==0 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_photocoag==1 & inclusion_died_prior==1

count if inclusion_mi_stroke==1 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_photocoag==1 & inclusion_died_prior==1

log close