/******************************************************************************* 

NAME: 			create_treatment_episodes.do

PURPOSE: 		calculate the duration of each drug in a sensible way 
				collapse overlapping drug prescriptions 
AUTHOR:			anna schultze
DATA USED: 		dataset prepared by Emma Powel in MSc/Yan-Ling Lu
DATA CREATED:   two stata datasets per drug, one in long and one in wide 

EDITS: 			n/a

*******************************************************************************/ 
/* HOUSEKEEPING===============================================================*/ 
global datadir "J:\EHR-Working\AnnaS\MSc\Yan-Ling_Lu\raw_data"
global tempdir "J:\EHR-Working\AnnaS\MSc\Yan-Ling_Lu\derived_data"
global cleandatadir "J:\EHR-Working\AnnaS\MSc\Yan-Ling_Lu\datafiles"

/* TREATMENT EPISODES=========================================================*/ 

foreach drug in Warfarin Apixaban Rivaroxaban Edoxaban Dabigatran { 
	
	di "`drug'"

	use $tempdir/treatment_episodes_temp_`drug', clear

	* create treatment episodes using the make_episodes program 
	make_episodes, patid(patid) start_date(rxst) end_date(rxend) allowable_gap(56) overlap(ignore) 
	
	* basic treatment episode descriptions 
	bysort patid: gen total_episodes = _N 
	bysort patid: egen first_rxend = min(end)
	format first_rxend %td 
	
	gen episode_duration = (end - start)
	bysort patid(start): gen time_between_eps = start[_n] - end[_n-1]
	
	* describe data 
	
	summ total_episodes, d
	summ episode_duration, d
	summ episode_duration if end == first_rxend, d 
	summ time_between_eps, d

	* clean up data
	sort patid start 
	keep patid start end actual_end total_episodes first_rxend
	rename end rxend
	rename start rxstart 
	rename actual_end actual_rxend
	
	* save as long 
	save $tempdir/treatment_episodes_long_`drug', replace
	
	* extract first discontinuation and save as 1 row pp 
	preserve
	count
	bysort patid(rxstart): keep if _n == 1 
	count
	save $tempdir/treatment_episodes_`drug', replace
	restore 
	
}