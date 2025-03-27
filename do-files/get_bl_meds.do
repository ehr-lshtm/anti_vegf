** Program to identify medications at baseline

* Author: Ruth Costello 
* Date: 11/11/2024

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
log using "$projdir\logs\2_meds_$S_DATE.log", append

run  "$dodir\pr_make_episodes.do"

* Baseline medications for those included in the study
foreach grp in antivegf cataract photocoag {
	
	display in red "*******************Drug issue file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\drugissue_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
	
	* merge diabetics drug codelist as different names and multiple drugs in one codelist
	merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 
	
	foreach script in acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin nsaid  {
	
		merge m:1 prodcodeid using "$savedir\codelists\med_`script'_codes_aurum", keepusing(drug_`script') nogen 
	
	}
	egen keep = rowmax(sub_mtf-drug_nsaid)
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
		
		merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
		
		* merge diabetics drug codelist as different names and multiple drugs in one codelist
	merge m:1 prodcodeid using "$savedir\codelists\med_antidiabetics_grouped_codes_aurum", keepusing(sub_*) nogen 
		
		foreach script in acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin nsaid {
	
		merge m:1 prodcodeid using "$savedir\codelists\med_`script'_codes_aurum", nogen keepusing(drug_`script')
	
		}
		egen keep = rowmax(sub_acarbose-drug_nsaid)
		tab keep 
		keep if keep!=.

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}
	compress 
	save "$savedir\\`grp'\all_drugs", replace
}

* Next need to identify which drugs people are taking at index. To do this will create treatment episodes from 6 months prior to index to 6 months after index. Prescriptions after index could be helpful to determine if prescriptions are ongoing.
* First generate datasets for each drug and clean duration variable
foreach grp in antivegf cataract photocoag {
	foreach i in  14 30 60 {
		use "$savedir\\`grp'\all_drugs", clear
		drop if patid==""
		merge m:1 dosageid using "$rawdata\common_dosages"
		drop if _merge==2
		drop _merge keep
		gen start_date = date(issuedate, "DMY")
		* Drop missing event dates
		drop if start_date==.
		* Keep records in 6 months prior to index and 6 months after index 
		gen days_prior_index = start_date-index_date 
		gen presc_keep = (days_prior_index>-183 & days_prior_index<=183)
		tab presc_keep 
		keep if presc_keep 
		destring daily_dose, replace
		* To determine duration of prescription: calculate duration based on quantity and daily_dose, drop calculated duration if 1 or less. Update calcuated duration to include duration from duration variable if missing when calculated - i.e. quantity and/or daily_dose is missing. 
		* Diabetes drugs
		foreach drug in acarbose first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i {
			preserve 
			tab sub_`drug' if sub_`drug'==1
			if r(N)>0 {
				keep if sub_`drug'==1
				sum quantity daily_dose duration, d
				gen duration_calc = quantity/daily_dose 
				replace duration_calc=. if duration_calc<=1
				replace duration_calc=duration if duration_calc==. & duration!=1 & duration!=0
				sum duration_calc duration 
				* Across all drugs the 99% duration frequently 92 days, one was 140 days. Assume prescriptions up to 4 months, therefore set durations over 4 months (120 days) to missing  
				replace duration_calc=. if duration_calc>120
				sum duration_calc
				di "Number of records"
				count 
				* Set missing durations to 28 days as this is median generally
				replace duration_calc=28 if duration_calc==.
				gen end_date = start_date+duration_calc  
				rename sub_`drug' drug_`drug'
				keep patid start_date end_date drug_`drug'
				* Create episodes of treatment
				make_episodes, patid(patid) start_date(start_date) end_date(end_date) allowable_gap(`i') overlap(ignore)
				* Merge on index date to determine if on treatment at index 
				merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid index_date) keep(match) nogen
				gen index_`drug' = (index_date>=start & index_date<=end)
				tab index_`drug' if index_`drug'==1
				if r(N)>0 {
					keep if index_`drug'
					keep patid index_`drug'
					duplicates drop 
					save "$savedir\\`grp'\cr_`drug'_index_`i'_prep", replace 
					restore		
				}
				else if r(N)==0 {
					restore 
				}
			}
			else if r(N)==0 {
				restore
			}
			}
		* Other drugs
		foreach drug in acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin nsaid {	
			preserve 
			tab drug_`drug' if drug_`drug'==1
			if r(N)>0 {
				keep if drug_`drug'==1
				sum quantity daily_dose duration, d
				gen duration_calc = quantity/daily_dose 
				replace duration_calc=. if duration_calc<=1
				replace duration_calc=duration if duration_calc==. & duration!=1 & duration!=0
				sum duration_calc duration 
				* Across all drugs the 99% duration frequently 92 days, one was 140 days. Assume prescriptions up to 4 months, therefore set durations over 4 months (120 days) to missing  
				replace duration_calc=. if duration_calc>120
				sum duration_calc
				di "Number of records"
				count 
				* Set missing durations to 28 days as this is median generally
				replace duration_calc=28 if duration_calc==.
				gen end_date = start_date+duration 
				keep patid start_date end_date drug_`drug'
				* Create episodes of treatment
				make_episodes, patid(patid) start_date(start_date) end_date(end_date) allowable_gap(`i') overlap(ignore)
				* Merge on index date to determine if on treatment at index 
				merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid index_date) keep(match) nogen
				gen index_`drug' = (index_date>=start & index_date<=end)
				tab index_`drug' if index_`drug'==1
				if r(N)>0 {
					keep if index_`drug'
					keep patid index_`drug'
					duplicates drop 
					save "$savedir\\`grp'\cr_`drug'_index_`i'_prep", replace 
					restore
				}
				else if r(N)==0 {
					restore 
				}
			}
			else if r(N)==0 {
				restore
			}
		}	
	}
}

foreach grp in antivegf cataract photocoag {
	foreach i in 14 30 60 {
		di `i'
		use "$savedir\\`grp'\cr_acarbose_index_`i'_prep", clear
		foreach drug in first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin nsaid {
			capture confirm file "$savedir\\`grp'\cr_`drug'_index_`i'_prep.dta"
			if   c(rc) {
										  di "`drug' file does not exist. Skipping..."
			}
			else {
										  di "Merging"
			merge 1:1 patid using "$savedir\\`grp'\cr_`drug'_index_`i'_prep", nogen 
			}
		}
		count 
		merge 1:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid) nogen
		foreach drug in acarbose first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin nsaid {
			capture confirm variable index_`drug'
			if !_rc {
				replace index_`drug'=0 if index_`drug'==.
			}
		}
		save "$savedir\\`grp'\cr_drugs_`i'_index", replace
	}
}


/*	OLD CODE	
	keep patid drug_*
	duplicates drop 
	codebook patid 
	count
	merge 1:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid) nogen
	codebook patid 
	foreach drug in acarbose first_gen_su second_gen_su third_gen_su dpp4i sglt2i glp1 glp2 guar mtf tzd glin ins mtf_dpp4i mtf_tzd mtf_sglt2i ins_glp1 3su_tzd sglt2i_dpp4i acei antiplatelet arb arni betablocker ccb loop_diuretic mra oac otherantihypertensive ppi statin {
		replace drug_`drug'=0 if drug_`drug'==.
	}
	save "$savedir\\`grp'\cr_bl_meds", replace
	
}
