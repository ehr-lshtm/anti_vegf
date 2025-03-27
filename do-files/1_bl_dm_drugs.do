** Program to identify diabetes medications at baseline for inclusion/exclusion

* Author: Ruth Costello 
* Date: 20/11/2024

/* 
What the program does:
- for each group and condition:
- Loops through each observation file 
- Merge on file of patients included in the study population
- Merge on medications codelist and keep matches
- Keep records prior to index
- generate flag for having medication at index if a record in 6 months prior to index 
- merge back on patient list and update flag to 0 if no records 
- generates file for each condition, then merges all together to create file of all conditions
*/

*********************************
log using "$projdir\logs\1_dm_meds_$S_DATE.log", append

* Baseline diabetes medications for inclusion/exclusion criteria 
foreach grp in cataract photocoag {
	
	display in red "*******************Drug issue file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\drugissue_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\prelim_pts_all", keepusing(patid yob first_code) keep(match) nogen
	
	* merge diabetics drug codelist as different names and multiple drugs in one codelist
	merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 
	
	egen keep = rowmax(sub_mtf-sub_sglt2i_dpp4i)
	tab keep 
	keep if keep!=.
	
	save `tempfile', replace
	
		/*******************************************************************************
	#A2. Loop through subsequent (from 2 onwards) separate test extract files in 
		turn and append the results to the first extract file saved in #1
	*******************************************************************************/
	if "`grp'"=="antivegf" {
		local j 5
	}
	if "`grp'"=="cataract" {
		local j 72
	}
	if "`grp'"=="photocoag" {
		local j 3
	}
	di `j'
	forvalues n=2/`j' {
		display in red "*******************Drug issue file number: `n'*******************"

		use "$rawdata\\`grp'\drugissue_`n'", clear
		
		merge m:1 patid using "$savedir\\`grp'\prelim_pts_all", keepusing(patid yob first_code) keep(match) nogen
		
		* merge diabetics drug codelist as different names and multiple drugs in one codelist
	merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 
		
		egen keep = rowmax(sub_mtf-sub_sglt2i_dpp4i)
		tab keep 
		keep if keep!=.

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}
	
	gen eventdate = date(issuedate, "DMY")
	* Drop missing event dates
	drop if eventdate==.
	* Keep records in 6 months prior to index
	gen days_prior_index = eventdate-first_code 
	gen prior_index = (days_prior_index>-183 & days_prior_index<=0)
	tab prior_index 
	keep if prior_index
	* Update so indicator on each row 
	* Diabetes drugs
	foreach drug in acarbose first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 guar mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i {
		bys patid: egen max_`drug' = max(sub_`drug')
		replace max_`drug' = 0 if max_`drug'==.
		tab max_`drug' sub_`drug'
		drop sub_`drug'
		rename max_`drug' drug_`drug'
	}
	
	keep patid drug_*
	duplicates drop 
	codebook patid 
	count
	
	save "$savedir\\`grp'\cr_dm_meds_prelim", replace
}

* Anti-vegf cohort

tempfile tempfile 
use "$rawdata\antivegf\drugissue_1", clear

merge m:1 patid using "$savedir\antivegf\prelim_pts", keepusing(patid yob first_code) keep(match) nogen

* merge diabetics drug codelist as different names and multiple drugs in one codelist
merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 

egen keep = rowmax(sub_mtf-sub_sglt2i_dpp4i)
tab keep 
keep if keep!=.

save `tempfile', replace

	/*******************************************************************************
#A2. Loop through subsequent (from 2 onwards) separate test extract files in 
	turn and append the results to the first extract file saved in #1
*******************************************************************************/
forvalues n=2/5 {
	display in red "*******************Drug issue file number: `n'*******************"

	use "$rawdata\antivegf\drugissue_`n'", clear
	
	merge m:1 patid using "$savedir\antivegf\prelim_pts", keepusing(patid yob first_code) keep(match) nogen
	
	* merge diabetics drug codelist as different names and multiple drugs in one codelist
merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 
	
	egen keep = rowmax(sub_mtf-sub_sglt2i_dpp4i)
	tab keep 
	keep if keep!=.

	* add the file containing records for the specified comorbidity
	* to make one file containing all specified comorbidiy records for the
	* clinical extract specified
	append using `tempfile'
	
	* save
	save `tempfile', replace
}

gen eventdate = date(issuedate, "DMY")
* Drop missing event dates
drop if eventdate==.
* Keep records in 6 months prior to index
gen days_prior_index = eventdate-first_code 
gen prior_index = (days_prior_index>-183 & days_prior_index<=0)
tab prior_index 
keep if prior_index
* Update so indicator on each row 
* Diabetes drugs
foreach drug in acarbose first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 guar mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i {
	bys patid: egen max_`drug' = max(sub_`drug')
	replace max_`drug' = 0 if max_`drug'==.
	tab max_`drug' sub_`drug'
	drop sub_`drug'
	rename max_`drug' drug_`drug'
}

keep patid drug_*
duplicates drop 
codebook patid 
count

save "$savedir\antivegf\cr_dm_meds_prelim", replace