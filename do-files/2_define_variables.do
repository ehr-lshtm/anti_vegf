/*==============================================================================
DO FILE NAME:			2_define_variables.do
DATE: 					16/10/2024
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	defines characteristics at index 

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\2_define_vars_$S_DATE.log", append

* Keep only patients inlcuded in study
* Creates file with index date, yob and gender (and DM flag for antivegf)
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\prelim_pts_4", clear
	keep if inclusion_died_prior==1 & inclusion_mi_stroke==1 & inclusion_egfr==1 & inclusion_reg==1 & inclusion_accept==1 & inclusion_age==1 & inclusion_prior_rx_`grp'  
	drop inclusion* _merge
	rename first_code index_date 
	count 
	save "$savedir\\`grp'\cr_study_pop", replace
}

* Exposure - Frequency of anti-vegf injections 
use "$rawdata\linked_data_20250106\linked_hes_procedures_epi", clear
gen antivegf = (opcs == "C794" | opcs == "C893")
keep if antivegf 
merge m:1 patid using "$savedir\antivegf\cr_study_pop", keep(match) keepusing(index_date bl_dm)
codebook patid 
gen admitdate = date(admidate, "YMD")
drop admidate 
gen evdateA = date(evdate, "YMD")
* Checking first code vs index date provided by CPRD
bys patid: egen min_date = min(admitdate)
bys patid: egen min_dateA = min(evdateA)
count if min_dateA!=index_date
* Identifying last code 
bys patid: egen last_injection = max(evdateA)
* Check time between injections 
bys patid (evdateA): gen time_since_prior = evdateA - evdateA[_n-1] if patid==patid[_n-1]
sum time_since_prior, d
* Gap of 6 months
gen gap_6_injections_i = time_since_prior>183 & time_since_prior!=. 
bys patid: egen gap_6_injections = max(gap_6_injections_i)
bys patid: egen first_gap_date_i = min(evdateA) if gap_6_injections==1
bys patid: egen first_gap_date = max(first_gap_date_i)
* Determine number of injections per year 
gen yr_evdate = year(evdateA)
bys patid yr_evdate: egen number_antivegf_yr = total(antivegf)
* Determine since start of follow-up 
bys patid: egen first = min(year(evdateA))
gen yr_since_fu = yr_evdate - first 
keep patid index_date bl_dm yr_evdate yr_since_fu number_antivegf_yr last_injection gap_6_injections first_gap_date
duplicates drop 
save "$savedir\antivegf\injections", replace 

* Demographics
* Run external scripts for ethnicity, BMI and smoking
do "$dodir\pr_getethnicitystatus_Aurum.do"
do "$dodir\pr_getallbmirecords_Aurum.do"
run "$dodir\pr_getsmokingstatus_Aurum.do"

* Get data to define smoking 
* Import data
* Smoking codes from observation files
foreach grp in antivegf cataract photocoag {
	
	display in red "*******************Observation file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
	
	merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_smoking_codes", keep(match) keepusing(smokstatus) nogen
	
	
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
			
		merge m:1 medcodeid using "$savedir\codelists\cprd_aurum_smoking_codes", keep(match) keepusing(smokstatus) nogen

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}

	save "$savedir\\`grp'\cr_obs_smoking", replace
	
	* Nicotine replacement from drug issue files
	use "$rawdata\\`grp'\drugissue_1", clear
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
			
	merge m:1 prodcodeid using "$savedir\codelists\cprd_aurum_nicotine_codes", keep(match) keepusing(smokstatus) nogen
		
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
		display in red "*******************Observation file number: `n'*******************"

		use "$rawdata\\`grp'\drugissue_`n'", clear
		
		merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid yob index_date) keep(match) nogen
			
		merge m:1 prodcodeid using "$savedir\codelists\cprd_aurum_nicotine_codes", keep(match) keepusing(smokstatus) nogen

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}
	
	save "$savedir\\`grp'\cr_drugissue_smoking", replace
	
}

* Run smoking algorithm program
foreach grp in antivegf cataract photocoag {
    
 use "$savedir\\`grp'\cr_study_pop", clear
 keep patid index_date
 noi pr_getsmokingstatus_Aurum, obsfile("$savedir\\`grp'\cr_obs_smoking") ///
 therapyfile("$savedir\\`grp'\cr_drugissue_smoking") ///
 smokingstatusvar(smokstatus) index(index_date)
 save "$savedir\\`grp'\cr_smoke_update", replace
}

* Diabetes - have flag of type of diabetes through codes, but will also need medications to determine/corroborate type as some people have both t1 and t2 codes, and some have only unclear codes
foreach grp in antivegf photocoag {
	use "$savedir\\`grp'\dm", clear
	* Take out records prior to yob 
	gen obsdate1=date(obsdate,"DMY")
	format obsdate1 %dD/N/CY
	drop obsdate 
	rename obsdate1 obsdate
	drop if obsdate==.
	gen yr_obs = year(obsdate)
	count if yr_obs <= yob
	drop if yr_obs <= yob
	duplicates drop
	tab bl_dm, m
	* Identify diabetes prior to first injection 
	gen prior_first_code = obsdate<=first_code 
	tab prior_first_code 
	keep if prior_first_code
	* Identify first DM code and then time since first code at index_date
	bys patid: egen first_dm_code = min(obsdate)
	gen yrs_dm = (first_code - first_dm_code)/365.25
	sum yrs_dm, d 
	* Determine type of diabetes
	bys patid: egen max_dm_type = max(bl_dm_type)
	* Check if records have type 1 and type 2 diabetes recorded.
	gen dm_t1 = (bl_dm_type==1)
	bys patid: egen max_dm_t1 = max(dm_t1)
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(patid) 
	drop if _merge==1
	drop _merge 
	bys patid: gen n=_n
	tab max_dm_type max_dm_t1 if n==1
	replace bl_dm=0 if bl_dm==.
	* Merge on diabetes drugs to help determine diabetes type
	merge m:1 patid using "$savedir\\`grp'\cr_drugs_index", keep(match) nogen
	drop index_acei-index_statin
	order index_ins, before(index_acarbose) 
	* Determine number of T2 DM drugs (excluding insulin)
	if "`grp'" == "antivegf" {
		egen drug_t2dm = rowtotal(index_acarbose-index_mtf_sglt2i)
	}
	else if "`grp'" == "photocoag" {
		egen drug_t2dm = rowtotal(index_acarbose-index_sglt2i_dpp4i)
	}
	* Determine number of DM drugs (including insulin)
	if "`grp'" == "antivegf" {
		egen drug_dm_count = rowtotal(index_ins-index_mtf_sglt2i)
	}
	else if "`grp'" == "photocoag" {
		egen drug_dm_count = rowtotal(index_ins-index_sglt2i_dpp4i)
	}
	foreach comb in mtf_dpp4i mtf_tzd mtf_sglt2i sglt2i_dpp4i {
		capture confirm variable index_`comb' 
		if !_rc {
			tab index_`comb'
			replace drug_dm_count = drug_dm_count + 1 if index_`comb'==1
		}
	}
		
	* clinical input from Jen Lees to determine diabetes type: 
	* If have both T1 and T2 code and only insulin prescribed = T1
	gen dm_type_f = max_dm_type
	replace dm_type_f = 1 if (max_dm_type==2 & max_dm_t1==1 & index_ins==1 & drug_t2dm==0)
	* If have both T1 and T2 code and no insulin prescribed or no drugs at all prescribed = T2 and therefore remain 2
	* Use same rules if have no clear codes
	replace dm_type_f = 1 if (max_dm_type==0 & index_ins==1 & drug_t2dm==0)
	replace dm_type_f = 2 if (max_dm_type==0 & drug_t2dm>=1)
	tab dm_type_f max_dm_t1 if n==1 & bl_dm==1
	
	keep patid bl_dm dm_type_f yrs_dm drug_dm_count
	duplicates drop 
	count 
	codebook patid 
	save "$savedir\\`grp'\cr_bl_dm", replace 
}

* Generate closest egfr to index date
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\eGFR", clear
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	duplicates drop patid obsdate SCr unit, force
	* Check if there are multiple measures on the same day
	bys patid obsdate: gen total = _N
	tab total 
	bys patid obsdate: egen min_egfr = min(egfr)
	bys patid obsdate: egen max_egfr = max(egfr)
	gen diff_egfr = min_egfr - max_egfr
	
	gen prior_index = (obsdate <= index_date)
	keep if prior_index
	bys patid: egen closest_egfr_date = max(obsdate)
	* Check time between index and closest egfr 
	gen days = closest_egfr_date - index_date
	sum days
	* Check difference between multiple measures on the same day that are closest to index
	sum diff_egfr if obsdate==closest_egfr_date & total>1
	* Rules for mulitple measures on same day - biggest difference is 20 so take lowest measure
	keep if obsdate==closest_egfr_date
	bys patid: egen prior_egfr = min(egfr)
	keep patid prior_egfr
	duplicates drop 
	count 
	codebook patid 
	* categorise into ckd stages
	egen prior_egfr_cat= cut(prior_egfr), at(0, 15, 30, 45, 60, 5000)
	*label define EGFR 0"stage 5" 15"stage 4" 30"stage 3b" 45"stage 3a" 60"no CKD"
	label values prior_egfr_cat EGFR
	label var prior_egfr_cat "eGFR category calc without eth + DN fudge factor"
	
	* * recode with appropriate category as reference
	recode prior_egfr_cat 0=5 15=4 30=3 45=2 60=0, generate(prior_ckd)
	*label define ckd 0"no CKD" 2"stage 3a" 3"stage 3b" 4"stage 4" 5"stage 5"
	label values prior_ckd ckd
	label var prior_ckd "CKD stage calc without eth + DN fudge factor"
	
	save "$savedir\\`grp'\cr_bl_egfr", replace
}

* Get ACR values for all three groups 
do "$dodir/prog_getACR_Aurum.do"

* Generate closest ACR to index 
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\ACR", clear
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	duplicates drop patid obsdate ACR unit, force
	* Check if there are multiple measures on the same day
	bys patid obsdate: gen total = _N
	tab total 
	bys patid obsdate: egen min_acr = min(ACR)
	bys patid obsdate: egen max_acr = max(ACR)
	gen diff_acr = min_acr - max_acr
	
	gen prior_index = (obsdate <= index_date)
	keep if prior_index
	bys patid: egen closest_acr_date = max(obsdate)
	* Check time between index and closest acr 
	gen days = closest_acr_date - index_date
	sum days
	* Check difference between multiple measures on the same day that are closest to index
	sum diff_acr if obsdate==closest_acr_date & total>1
	* Rules for mulitple measures on same day - biggest difference is 20 so take lowest measure
	keep if obsdate==closest_acr_date
	bys patid: egen prior_acr = min(ACR)
	keep patid prior_acr
	duplicates drop 
	count 
	codebook patid 
	* categorise into ckd stages
	egen prior_acr_cat= cut(prior_acr), at(0, 30, 300, 5000) icodes
	label define ACR 0"A1 Normal to mildly increased" 1 "A2 Moderately increased" 2 "A3 Severely increased"
	label values prior_acr_cat ACR
	label var prior_acr_cat "Albuminuria category"
	
	save "$savedir\\`grp'\cr_bl_acr", replace
}

* Healthcare utilisation in year prior to index 
/* Referrals to secondary care - decided to park as Helen Strongman has investigated using this data but decided that it wasn't well coded and couldn't be determined very well
foreach grp in antivegf cataract photocoag {
	
	display in red "*******************Observation file number: 1*******************"
	tempfile tempfile 
	use "$rawdata\\`grp'\observation_1", clear
	
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	
	merge m:1 obsid using "$rawdata\referral_1", keep(match) nogen
	
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
		
		merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	
		merge m:1 obsid using "$rawdata\referral_1", keep(match) nogen


		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using `tempfile'
		
		* save
		save `tempfile', replace
	}
	
	* keep only those in year prior to index date 
	replace obsdate=enterdate if (obsdate=="" & enterdate!="")
	gen obsdate1=date(obsdate,"DMY")
	format obsdate1 %dD/N/CY

	drop obsdate 
	rename obsdate1 obsdate
	drop if obsdate==.
	*****************************/
	
	
* Number of creatinine measures in year prior 
foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\eGFR", clear 
	duplicates drop patid obsdate, force 
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", keepusing(index_date) keep(match) nogen
	codebook patid 
	gen days_prior = index_date - obsdate 
	gen yr_prior = days_prior>=0 & days_prior<=366
	keep if yr_prior 
	codebook patid 
	bys patid: gen scr_meas_yr_prior = _N
	sum scr_meas_yr_prior, d
	keep patid scr_meas_yr_prior
	duplicates drop 
	* Eligibility means creatinine measures could be more than one year prior to index, therefore some people could have no measures in last year
	merge m:1 patid using "$savedir\\`grp'\cr_study_pop", nogen
	replace scr_meas_yr_prior = 0 if scr_meas_yr_prior==.
	save "$savedir\\`grp'\cr_scr_meas_yr_prior", replace
}

* creating overall dataset
tempfile tempfile1 
use "$rawdata\antivegf\patient_1", clear
merge 1:1 patid using "$savedir\antivegf\cr_study_pop", keepusing(patid index_date) keep(match) nogen
gen group = 1
save `tempfile1'

tempfile tempfile2 
use "$rawdata\photocoag\patient_1", clear 
merge 1:1 patid using "$savedir\photocoag\cr_study_pop", keepusing(patid index_date) keep(match) nogen
gen group=2
save `tempfile2'

tempfile tempfile3 
use "$rawdata\cataract\patient_1", clear 
merge 1:1 patid using "$savedir\cataract\cr_study_pop", keepusing(patid index_date) keep(match) nogen
gen group=3
save `tempfile3'

use `tempfile1', clear 
append using `tempfile2'
append using `tempfile3'
tempfile tempfile4
save "$savedir\cr_all_study_pop", replace

* Number of ophthalmology outpatient appointments
use "$rawdata\linked_data_20250106\linked_hesop_appointment", clear 
merge 1:m patid attendkey using "$rawdata\linked_data_20250106\linked_hesop_clinical", nogen
destring mainspef tretspef, force gen(mainspefA tretspefA)
keep if (mainspefA==130 | mainspefA == 460 | tretspefA==130 | tretspefA== 460)
merge m:1 patid using "$savedir\cr_all_study_pop", keep(match) keepusing(index_date group) nogen
* Keeping only those appointments attended 
keep if (attended == "5" | attended == "6")
save "$savedir\cr_all_op_appt_attend", replace
gen apptdateA = date(apptdate, "YMD")
* Keep only year prior to index
gen days_prior = index_date - apptdateA 
gen yr_prior = days_prior>=0 & days_prior<=366
keep if yr_prior
bys patid: egen tot_appts_yr_prior = total(yr_prior)
sum tot_appts_yr_prior, d
keep patid tot_appts_yr_prior
duplicates drop 
codebook patid 
count
merge m:1 patid using "$savedir\cr_all_study_pop", keepusing(group) nogen
replace tot_appts_yr_prior=0 if tot_appts_yr_prior==.
preserve 
keep if group==1
count
save "$savedir\antivegf\cr_bl_op_appts"
restore 

preserve 
keep if group==2
count
save "$savedir\photocoag\cr_bl_op_appts"
restore

preserve 
keep if group==3
count
save "$savedir\cataract\cr_bl_op_appts"
restore

* Identify death dates from linked ONS data and compare to dates in CPRD
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
* Compare to CPRD death date
gen same_cprd = (cprd_death_date==date_of_death)
tab same_cprd
bys patid (date_of_death): gen diff_ons = date_of_death!=date_of_death[_n-1] & patid==patid[_n-1]
* Determine first death date
bys patid: egen death_date = min(date_of_death)
* Only two people who haev multiple ONS death dates and one of those dates matches CPRD death date - in both cases only one day out from minimum therefore keeping as minimum  
gen days = date_of_death-death_date
sum days if days!=0
keep patid death_date 
duplicates drop 
codebook patid 
count 
save "$savedir\linked_death_dates", replace 

do "$dodir\get_bl_meds.do"

* Explore impact of different windows on drugs at index 
foreach i in 14 30 60 {
	use "$savedir\antivegf\cr_drugs_`i'_index", clear
	merge 1:1 patid using "$savedir\antivegf\cr_bl_dm", gen(dm_merge) keep(match)
	keep if bl_dm==1
	gen antivegf=1
	append using "$savedir\photocoag\cr_drugs_`i'_index"
	replace antivegf = 0 if antivegf==.
	dtable i.index_*, by(antivegf) export($projdir\output\explore_window_`i'.xlsx, replace)
}

* Run script for comorbidities at index 
do "$dodir\pr_get_comorbidities_at_index.do"

foreach grp in antivegf cataract photocoag {
	use "$savedir\\`grp'\cr_study_pop", clear
	
	* Define end of follow-up (min of end of study period, death, end registration)
	merge 1:1 patid using "$rawdata\\`grp'\patient_1", keep(match) keepusing(regenddate cprd_ddate) nogen 
	merge 1:1 patid using "$savedir\linked_death_dates"
	drop if _merge==2
	drop _merge 
	gen reg_end = date(regenddate, "DMY")
	* Use CPRD death date while don't have ONS data
	gen cprd_death_date = date(cprd_ddate, "DMY")
	gen end_fu = min(reg_end, death_date, date("01Mar2020", "DMY"))
	format reg_end death_date cprd_death_date end_fu %dD/N/CY
	count if end_fu<index_date
	gen fu_time = end_fu - index_date
	sum fu_time, d
	
	** Demographics 
	* age at index 
	gen yr_index = year(index_date)
	gen age_at_index = yr_index - yob 
	sum age
	
	* Gender 
	label define sex 1 "Male" 2 "Female"
	label values gender sex
	
	* ethnicity 
	merge 1:1 patid using "$savedir\\`grp'\cr_ethnicity", gen(eth_merge) keep(match)
	
	* BMI 
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_bmi", gen(bmi_merge) 
	drop if bmi_merge==2
	* Heaviest ever man in UK had a BMI of 155 so take out any above 155
	replace bl_bmi = . if bl_bmi>155
	* Define BMI categories
	egen bmi_cat = cut(bl_bmi), at(0, 1, 18.5, 24.9, 29.9, 39.9, 100) icodes
	bys bmi_cat: sum bl_bmi
	* add missing . to zero category
	replace bmi_cat = 0 if bmi_cat==. 
	label define bmi 0 "Missing" 1 "Underweight" 2 "Healthy range" 3 "Overweight" 4 "Obese" 5 "Severely obese"
	label values bmi_cat bmi
	
	* Smoking 
	merge 1:1 patid using "$savedir\\`grp'\cr_smoke_update", gen(smok_merge) keepusing(smokstatus) 
	drop if smok_merge==2
	
	** Diabetes variables
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_dm", gen(dm_merge) keep(match)
	
	* Drugs at index
	merge 1:1 patid using "$savedir\\`grp'\cr_drugs_60_index", gen(meds_merge) keep(match)
	
	* Comorbidities at index
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_comorbidities_hes", gen(comorbid_merge) keep(match)
	
	* Generate duration of eye disease at index 
	gen yrs_retinopathy =  (index_date - first_code_retinopathy)/365.25
	replace yrs_retinopathy = 0 if yrs_retinopathy==.
	gen yrs_eye_dis = (index_date - first_code_eye_dis)/365.25
	replace yrs_eye_dis = 0 if yrs_eye_dis==.
	
	* eGFR at index 
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_egfr", gen(egfr_merge) keep(match)
	
	* Update kidney failure indicator to include egfr<15 
	replace bl_kidney_failure = 1 if prior_egfr < 15
	
	* ACR at index 
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_acr", gen(acr_merge) 
	drop if acr_merge==2
	
	* Ophthalmology appointments at index 
	merge 1:1 patid using "$savedir\\`grp'\cr_bl_op_appts", gen(op_merge) 
	
	* Creatinine measures at index 
	merge 1:1 patid using "$savedir\\`grp'\cr_scr_meas_yr_prior", gen(creatinine_merge)
	
	* IMD
	merge 1:1 patid using "$rawdata\linked_data_20250106\linked_imd", gen(imd_merge)
	drop if imd_merge==2
	destring e2019_imd_5, gen(imd)
	drop e2019_imd_5
	replace imd=0 if imd==.
	
	drop *_merge
	save "$savedir\\`grp'\cr_vars", replace
}


	