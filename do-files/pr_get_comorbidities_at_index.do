** Program to identify conditions at baseline

* Author: Ruth Costello 
* Date: 06/11/2024

/* 
What the program does:
- for each group and condition:
- Loops through each observation file 
- Merge on file of patients included in the study population
- Merge on condition codelist and keep matches
- Keep records prior to index
- generate flag for having conidition at index if any record prior to index 
- merge back on patient list and update flag to 0 if no records 
- generates file for each condition, then merges all together to create file of all conditions
*/

*********************************
log using "$projdir\logs\2_comorbidities_$S_DATE.log", append


foreach grp in antivegf cataract photocoag {
	
	display in red "*******************Observation file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
	
	foreach script in af amd depression hypertension kidney_failure mi ncm neph_syndrome retinopathy rrt rvo stroke zoster hf pad copd cancer neuropathy amputation {
	
		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_`script'", nogen keepusing(`script'_cprd)
	
	}
	egen keep = rowmax(af_cprd-amputation_cprd)
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
		
		merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
			
		foreach script in af amd depression hypertension kidney_failure mi ncm neph_syndrome retinopathy rrt rvo stroke zoster hf pad copd cancer neuropathy amputation  {
	
		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_codes_`script'", nogen keepusing(`script'_cprd)
	
		}
		egen keep = rowmax(af_cprd-amputation_cprd)
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
	* Keep records prior to index
	gen prior_index = eventdate<=index_date 
	tab prior_index 
	keep if prior_index
	
	* Check year of event isn't prior to yob 
	gen year_event = year(eventdate)
	sum year_event if year_event<yob
	drop if year_event<yob 
	
	* Spread vars so single row per person
	foreach script in af amd depression hypertension kidney_failure mi ncm neph_syndrome retinopathy rrt rvo stroke zoster hf pad copd cancer neuropathy amputation  {
			bys patid: egen bl_`script' = max(`script'_cprd)
			bys patid: egen first_code_`script' = min(eventdate) if `script'_cprd==1
			format first_code_`script' %dD/N/CY
			* Spread to all rows
			bys patid (first_code_`script'): replace first_code_`script'=first_code_`script'[1] if first_code_`script'==.
		}
	
	keep patid bl_af-first_code_amputation
	duplicates drop 
	codebook patid 
	merge 1:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid) nogen
	* Tidying 
	foreach script in af amd depression hypertension kidney_failure mi ncm neph_syndrome retinopathy rrt rvo stroke zoster hf pad copd cancer neuropathy amputation  {
		replace bl_`script'=0 if bl_`script'==.
	}
	codebook patid 
	save "$savedir\\`grp'\cr_bl_comorbidities", replace
	
}

* Comorbidities from HES data

use "$rawdata\linked_data_20250106\linked_hes_diagnosis_hosp", clear 
drop if icd=="-."
* Eye disease
* Retinopathy
merge m:1 icd using "$savedir\codelists\hes_eye_dis_codes", nogen
destring hes_eye_dis, replace
gen hes_retinopathy = (hes_eye_dis==2)

* Amputation
gen hes_amputation = icd == "Y83.5"

* Nephrotic syndrome
gen hes_neph_syndrome = icd == "N04.0"

* MI 
merge m:1 icd using "$savedir\codelists\hes_mi_codes", nogen
drop if patid==""

* Stroke 
merge m:1 icd using "$savedir\codelists\hes_stroke_codes", nogen
drop if patid==""

* HF
merge m:1 icd using "$savedir\codelists\hes_hf", nogen
drop if patid==""
rename hf_hes hes_hf

* Arrhythmia
merge m:1 icd using "$savedir\codelists\hes_arrhythmia_codes", nogen
drop if patid==""

* Cancer
* Cancer doesn't use decimal in codes - take first 3 characters from icd variable
gen icd_3 = substr(icd, 1,3)
merge m:1 icd_3 using "$savedir\codelists\hes_cancer_codes", nogen
drop if patid==""

* Depression
merge m:1 icd_3 using "$savedir\codelists\hes_depression", nogen keepusing(hes_depression)
drop if patid==""
		
* COPD
merge m:1 icd using "$savedir\codelists\hes_copd_codes", nogen
drop if patid==""		

* Kidney failure 
merge m:1 icd using "$savedir\codelists\hes_kidney_failure_codes", nogen
drop if patid==""	

* Zoster
merge m:1 icd using "$savedir\codelists\hes_zoster", nogen
drop if patid==""
rename zoster_hes hes_zoster

foreach var in eye_dis amputation mi stroke hf arrhythmia cancer depression copd kidney_failure zoster {
	replace hes_`var'=0 if hes_`var'==.
	tab hes_`var'
}
* Keep only records where there are records for any disease 
drop icd_description icd_3
egen hes_any = rowmax(hes_eye_dis-hes_zoster)
foreach var in eye_dis amputation neph_syndrome mi stroke hf arrhythmia cancer depression copd kidney_failure zoster {
	count if hes_`var'!=0 & hes_any==0
}
drop if hes_any==0
gen admitdate = date(admidate, "YMD")
drop admidate 
save "$savedir/all_hosp_bl", replace
* Create file for each cohort
foreach grp in antivegf cataract photocoag {
	preserve 
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	* Update eye disease variable as have another variable for retinopathy 
	replace hes_eye_dis =0 if hes_eye_dis==2
	foreach var in eye_dis retinopathy amputation neph_syndrome mi stroke hf arrhythmia cancer depression copd kidney_failure {
		* Determine if row is prior to index 
		gen prior_`var'_i = (admitdate<=index_date & hes_`var'==1)
		* Drop records after index 
		drop if prior_`var'_i==0 & hes_`var'!=0
		* Create static flag for prior comorbidities
		bys patid: egen prior_`var' = max(prior_`var'_i)
		* Identify first code for each comorbidity
		bys patid: egen first_code_`var' = min(admitdate) if prior_`var'_i==1
		bys patid: egen first_code_`var'_hes = max(first_code_`var')
		drop first_code_`var'
	}
	count if admitdate>index_date
	keep patid prior_eye_dis prior_retinopathy prior_amputation prior_neph_syndrome prior_mi prior_stroke prior_hf prior_arrhythmia prior_cancer prior_depression prior_copd prior_kidney_failure first_code*
	duplicates drop 
	count 
	codebook patid 
	save "$savedir\\`grp'\cr_hes_bl", replace
	restore 
}

* Merge together with comorbidities in GP records
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\cr_bl_comorbidities", clear
	merge 1:1 patid using "$savedir\\`grp'\cr_hes_bl"
	* Edit some variable names to match Aurum variables 
	gen bl_eye_dis = (bl_amd==1 | bl_rvo==1)
	egen first_code_eye_dis = rowmin(first_code_amd first_code_rvo)
	rename prior_arrhythmia prior_af 
	rename first_code_arrhythmia_hes first_code_af_hes
	foreach var in eye_dis retinopathy amputation neph_syndrome mi stroke hf af cancer depression copd kidney_failure {
		replace bl_`var' = prior_`var' if bl_`var'==0 & prior_`var'==1
		replace first_code_`var' = first_code_`var'_hes if first_code_`var'_hes<first_code_`var'
		drop prior_`var' first_code_`var'_hes 
	}
	save "$savedir\\`grp'\cr_bl_comorbidities_hes", replace
}	
		
	