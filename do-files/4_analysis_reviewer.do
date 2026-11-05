/*==============================================================================
DO FILE NAME:			4_analysis_reviewer.do
DATE: 					22/10/2026
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	new cumulative incidence method

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020
This script refines the study populations based on the inclusion/exclusion criteria.
==============================================================================*/

log using "$projdir\logs\4_analysis_reviewer_$S_DATE.log", append

* Using stpm3 and estimating cummulative incidence taking into account competing risk of death
* Diabetes group

use "$savedir\an_dm_main_analysis_ps", clear
file open tablech using "$projdir/output/cumulative_hazard_cr.txt", write text replace
file write tablech ("Group") _tab ("Outcome") _tab ("Exposure") _tab ("Cumulative incidence") _tab ("95% confidence interval") _tab _tab ("Difference") _tab ("95% Confidence interval") _n

local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'

forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'

	* Create competing risks indicator: 1=outcome, 2=death, 0=censored
	cap drop `a'_cr
	gen `a'_cr = `a'
	replace `a'_cr = 2 if death_date == `b' & `a' == 0
	
	* Model 1: Primary outcome cause-specific hazard
	stset `b' [pweight=att_weight], failure(`a'_cr==1) origin(index_date) enter(index_date) id(patid) scale(365.25)
	stpm3 antivegf, df(2) scale(lncumhazard) eform
	estimates store outcome

	* Model 2: Death cause-specific hazard
	stset `b' [pweight=att_weight], failure(`a'_cr==2) origin(index_date) enter(index_date) id(patid) scale(365.25)
	stpm3 antivegf, df(2) scale(lncumhazard) eform
	estimates store death

	* Time variable
	summ _t
	local tmax = r(max)
	cap drop years
	range years 0 `tmax' 100

	* Standardised cause-specific CIFs
	* at1=photocoagulation (reference), at2=antivegf, so cif_diff = antivegf - photocoagulation
	standsurv, crmodels(outcome death) cif ci timevar(years) frame(cif2, replace) ///
		at1(antivegf 0) at2(antivegf 1) atvar(F_laser F_antivegf) ///
		contrast(difference) contrastvar(cif_diff)

	* Scale to percentages and extract end-of-follow-up values within frame
	frame cif2 {
		foreach v of varlist F_laser_outcome F_laser_outcome_lci F_laser_outcome_uci ///
		                     F_antivegf_outcome F_antivegf_outcome_lci F_antivegf_outcome_uci ///
		                     cif_diff_outcome cif_diff_outcome_lci cif_diff_outcome_uci {
			replace `v' = 100 * `v'
		}

		sum years
		local new_max = r(max)

		sum F_laser_outcome if years == `new_max'
		local cif_la = r(mean)
		sum F_laser_outcome_lci if years == `new_max'
		local cif_la_lci = r(mean)
		sum F_laser_outcome_uci if years == `new_max'
		local cif_la_uci = r(mean)

		sum F_antivegf_outcome if years == `new_max'
		local cif_av = r(mean)
		sum F_antivegf_outcome_lci if years == `new_max'
		local cif_av_lci = r(mean)
		sum F_antivegf_outcome_uci if years == `new_max'
		local cif_av_uci = r(mean)

		sum cif_diff_outcome if years == `new_max'
		local cif_diff_val = r(mean)
		sum cif_diff_outcome_lci if years == `new_max'
		local cif_diff_lci = r(mean)
		sum cif_diff_outcome_uci if years == `new_max'
		local cif_diff_uci = r(mean)
	}

	* Write photocoagulation row
	file write tablech ("Diabetes") _tab ("`a'") _tab ("Photocoagulation") _tab ///
		(`cif_la') _tab (`cif_la_lci') _tab (`cif_la_uci') _tab _n

	* Write antivegf row with difference
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab ///
		(`cif_av') _tab (`cif_av_lci') _tab (`cif_av_uci') _tab ///
		(`cif_diff_val') _tab (`cif_diff_lci') _tab (`cif_diff_uci') _tab _n

	* Plot
	frame cif2 {
		twoway  (rarea F_laser_outcome_lci F_laser_outcome_uci years, color(red%25)) ///
				(rarea F_antivegf_outcome_lci F_antivegf_outcome_uci years, color(blue%25)) ///
				(line F_laser_outcome years, sort lcolor(red)) ///
				(line F_antivegf_outcome years, sort lcolor(blue) lpattern(dash)) ///
				, legend(order(1 "Photocoagulation" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
				title("A", justification(left) size(med)) ///
				yscale(range(0, 35)) ///
				ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small)) ///
				xlabel(0 (1) 7.3, labsize(small)) ///
				ytitle("Cause-specific CIF (%)", size(medsmall)) ///
				xtitle("Years since index date", size(medsmall)) ///
				graphregion(fcolor(white)) name(dm_`a', replace)
	}

	graph export "$projdir/output/adjcurv_dm_`a'_cr.tif", as(tif) width(1271) height(762) replace

	* Cleanup
	cap drop `a'_cr years
	frame drop cif2
	estimates drop outcome death
}

drop if bl_hypertension==1

* Create competing risks indicator: 1=outcome, 2=death, 0=censored
cap drop event_hypertension_cr
gen event_hypertension_cr = event_hypertension
replace event_hypertension_cr = 2 if death_date == end_hypertension & event_hypertension == 0

* Model 1: Primary outcome cause-specific hazard
stset end_hypertension [pweight=att_weight], failure(event_hypertension_cr==1) origin(index_date) enter(index_date) id(patid) scale(365.25)
stpm3 antivegf, df(2) scale(lncumhazard) eform
estimates store outcome

* Model 2: Death cause-specific hazard
stset end_hypertension [pweight=att_weight], failure(event_hypertension_cr==2) origin(index_date) enter(index_date) id(patid) scale(365.25)
stpm3 antivegf, df(2) scale(lncumhazard) eform
estimates store death

* Time variable
summ _t
local tmax = r(max)
cap drop years
range years 0 `tmax' 100

* Standardised cause-specific CIFs
* at1=photocoagulation (reference), at2=antivegf, so cif_diff = antivegf - photocoagulation 
standsurv, crmodels(outcome death) cif ci timevar(years) frame(cif2, replace) ///
	at1(antivegf 0) at2(antivegf 1) atvar(F_laser F_antivegf) ///
	contrast(difference) contrastvar(cif_diff)

* Scale to percentages and extract end-of-follow-up values within frame
frame cif2 {
	foreach v of varlist F_laser_outcome F_laser_outcome_lci F_laser_outcome_uci ///
						 F_antivegf_outcome F_antivegf_outcome_lci F_antivegf_outcome_uci ///
						 cif_diff_outcome cif_diff_outcome_lci cif_diff_outcome_uci {
		replace `v' = 100 * `v'
	}

	sum years
	local new_max = r(max)

	sum F_laser_outcome if years == `new_max'
	local cif_la = r(mean)
	sum F_laser_outcome_lci if years == `new_max'
	local cif_la_lci = r(mean)
	sum F_laser_outcome_uci if years == `new_max'
	local cif_la_uci = r(mean)

	sum F_antivegf_outcome if years == `new_max'
	local cif_av = r(mean)
	sum F_antivegf_outcome_lci if years == `new_max'
	local cif_av_lci = r(mean)
	sum F_antivegf_outcome_uci if years == `new_max'
	local cif_av_uci = r(mean)

	sum cif_diff_outcome if years == `new_max'
	local cif_diff_val = r(mean)
	sum cif_diff_outcome_lci if years == `new_max'
	local cif_diff_lci = r(mean)
	sum cif_diff_outcome_uci if years == `new_max'
	local cif_diff_uci = r(mean)
}

* Write photocoagulation row
file write tablech ("Diabetes") _tab ("event_hypertension") _tab ("Photocoagulation") _tab ///
	(`cif_la') _tab (`cif_la_lci') _tab (`cif_la_uci') _tab _n

* Write antivegf row with difference
file write tablech _tab ("event_hypertension") _tab ("Antivegf") _tab ///
	(`cif_av') _tab (`cif_av_lci') _tab (`cif_av_uci') _tab ///
	(`cif_diff_val') _tab (`cif_diff_lci') _tab (`cif_diff_uci') _tab _n

* Plot
frame cif2 {
	twoway  (rarea F_laser_outcome_lci F_laser_outcome_uci years, color(red%25)) ///
			(rarea F_antivegf_outcome_lci F_antivegf_outcome_uci years, color(blue%25)) ///
			(line F_laser_outcome years, sort lcolor(red)) ///
			(line F_antivegf_outcome years, sort lcolor(blue) lpattern(dash)) ///
			, legend(order(1 "Photocoagulation" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
			title("C", justification(left) size(med)) ///
			yscale(range(0, 35)) ///
			ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small)) ///
			xlabel(0 (1) 7.3, labsize(small)) ///
			ytitle("Cause-specific CIF (%)", size(medsmall)) ///
			xtitle("Years since index date", size(medsmall)) ///
			graphregion(fcolor(white)) name(dm_event_hypertension, replace)
}

graph export "$projdir/output/adjcurv_dm_event_hypertension_cr.tif", as(tif) width(1271) height(762) replace

* No diabetes cohort 

use "$savedir\an_nodm_main_analysis_ps", clear

local event "egfr_40 egfr_40_sustained acr_increased event_af event_kidney_failure_15 event_mi event_neph_syndrome event_stroke event_hf event_pad event_cv event_cvd_death event_zoster"
local end_date "end_egfr_40 end_egfr_40_sustained end_acr end_af end_kidney_failure_15 end_mi end_neph_syndrome end_stroke end_hf end_pad end_cv end_cvd_death end_zoster"
local n: word count `event'

forvalues i=1/`n' {
	local a : word `i' of `event'
	local b : word `i' of `end_date'

	* Create competing risks indicator: 1=outcome, 2=death, 0=censored
	cap drop `a'_cr
	gen `a'_cr = `a'
	replace `a'_cr = 2 if death_date == `b' & `a' == 0
	
	* Model 1: Primary outcome cause-specific hazard
	stset `b' [pweight=att_weight], failure(`a'_cr==1) origin(index_date) enter(index_date) id(patid) scale(365.25)
	stpm3 antivegf, df(2) scale(lncumhazard) eform
	estimates store outcome

	* Model 2: Death cause-specific hazard
	stset `b' [pweight=att_weight], failure(`a'_cr==2) origin(index_date) enter(index_date) id(patid) scale(365.25)
	stpm3 antivegf, df(2) scale(lncumhazard) eform
	estimates store death

	* Time variable
	summ _t
	local tmax = r(max)
	cap drop years
	range years 0 `tmax' 100

	* Standardised cause-specific CIFs
	* at1=photocoagulation (reference), at2=antivegf, so cif_diff = antivegf - photocoagulation
	standsurv, crmodels(outcome death) cif ci timevar(years) frame(cif2, replace) ///
		at1(antivegf 0) at2(antivegf 1) atvar(F_laser F_antivegf) ///
		contrast(difference) contrastvar(cif_diff)

	* Scale to percentages and extract end-of-follow-up values within frame
	frame cif2 {
		foreach v of varlist F_laser_outcome F_laser_outcome_lci F_laser_outcome_uci ///
		                     F_antivegf_outcome F_antivegf_outcome_lci F_antivegf_outcome_uci ///
		                     cif_diff_outcome cif_diff_outcome_lci cif_diff_outcome_uci {
			replace `v' = 100 * `v'
		}

		sum years
		local new_max = r(max)

		sum F_laser_outcome if years == `new_max'
		local cif_la = r(mean)
		sum F_laser_outcome_lci if years == `new_max'
		local cif_la_lci = r(mean)
		sum F_laser_outcome_uci if years == `new_max'
		local cif_la_uci = r(mean)

		sum F_antivegf_outcome if years == `new_max'
		local cif_av = r(mean)
		sum F_antivegf_outcome_lci if years == `new_max'
		local cif_av_lci = r(mean)
		sum F_antivegf_outcome_uci if years == `new_max'
		local cif_av_uci = r(mean)

		sum cif_diff_outcome if years == `new_max'
		local cif_diff_val = r(mean)
		sum cif_diff_outcome_lci if years == `new_max'
		local cif_diff_lci = r(mean)
		sum cif_diff_outcome_uci if years == `new_max'
		local cif_diff_uci = r(mean)
	}

	* Write photocoagulation row
	file write tablech ("No diabetes") _tab ("`a'") _tab ("Cataract") _tab ///
		(`cif_la') _tab (`cif_la_lci') _tab (`cif_la_uci') _tab _n

	* Write antivegf row with difference
	file write tablech _tab ("`a'") _tab ("Antivegf") _tab ///
		(`cif_av') _tab (`cif_av_lci') _tab (`cif_av_uci') _tab ///
		(`cif_diff_val') _tab (`cif_diff_lci') _tab (`cif_diff_uci') _tab _n

	* Plot
	frame cif2 {
		twoway  (rarea F_laser_outcome_lci F_laser_outcome_uci years, color(red%25)) ///
				(rarea F_antivegf_outcome_lci F_antivegf_outcome_uci years, color(blue%25)) ///
				(line F_laser_outcome years, sort lcolor(red)) ///
				(line F_antivegf_outcome years, sort lcolor(blue) lpattern(dash)) ///
				, legend(order(1 "Cataract" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
				title("B", justification(left) size(med)) ///
				yscale(range(0, 35)) ///
				ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small)) ///
				xlabel(0 (1) 7.3, labsize(small)) ///
				ytitle("Cause-specific CIF (%)", size(medsmall)) ///
				xtitle("Years since index date", size(medsmall)) ///
				graphregion(fcolor(white)) name(nodm_`a', replace)
	}

	graph export "$projdir/output/adjcurv_nodm_`a'_cr.tif", as(tif) width(1271) height(762) replace

	* Cleanup
	cap drop `a'_cr years
	frame drop cif2
	estimates drop outcome death
}

drop if bl_hypertension==1

* Create competing risks indicator: 1=outcome, 2=death, 0=censored
cap drop event_hypertension_cr
gen event_hypertension_cr = event_hypertension
replace event_hypertension_cr = 2 if death_date == end_hypertension & event_hypertension == 0

* Model 1: Primary outcome cause-specific hazard
stset end_hypertension [pweight=att_weight], failure(event_hypertension_cr==1) origin(index_date) enter(index_date) id(patid) scale(365.25)
stpm3 antivegf, df(2) scale(lncumhazard) eform
estimates store outcome

* Model 2: Death cause-specific hazard
stset end_hypertension [pweight=att_weight], failure(event_hypertension_cr==2) origin(index_date) enter(index_date) id(patid) scale(365.25)
stpm3 antivegf, df(2) scale(lncumhazard) eform
estimates store death

* Time variable
summ _t
local tmax = r(max)
cap drop years
range years 0 `tmax' 100

* Standardised cause-specific CIFs
* at1=photocoagulation (reference), at2=antivegf, so cif_diff = antivegf - photocoagulation 
standsurv, crmodels(outcome death) cif ci timevar(years) frame(cif2, replace) ///
	at1(antivegf 0) at2(antivegf 1) atvar(F_laser F_antivegf) ///
	contrast(difference) contrastvar(cif_diff)

* Scale to percentages and extract end-of-follow-up values within frame
frame cif2 {
	foreach v of varlist F_laser_outcome F_laser_outcome_lci F_laser_outcome_uci ///
						 F_antivegf_outcome F_antivegf_outcome_lci F_antivegf_outcome_uci ///
						 cif_diff_outcome cif_diff_outcome_lci cif_diff_outcome_uci {
		replace `v' = 100 * `v'
	}

	sum years
	local new_max = r(max)

	sum F_laser_outcome if years == `new_max'
	local cif_la = r(mean)
	sum F_laser_outcome_lci if years == `new_max'
	local cif_la_lci = r(mean)
	sum F_laser_outcome_uci if years == `new_max'
	local cif_la_uci = r(mean)

	sum F_antivegf_outcome if years == `new_max'
	local cif_av = r(mean)
	sum F_antivegf_outcome_lci if years == `new_max'
	local cif_av_lci = r(mean)
	sum F_antivegf_outcome_uci if years == `new_max'
	local cif_av_uci = r(mean)

	sum cif_diff_outcome if years == `new_max'
	local cif_diff_val = r(mean)
	sum cif_diff_outcome_lci if years == `new_max'
	local cif_diff_lci = r(mean)
	sum cif_diff_outcome_uci if years == `new_max'
	local cif_diff_uci = r(mean)
}

* Write photocoagulation row
file write tablech ("No diabetes") _tab ("event_hypertension") _tab ("Cataract") _tab ///
	(`cif_la') _tab (`cif_la_lci') _tab (`cif_la_uci') _tab _n

* Write antivegf row with difference
file write tablech _tab ("event_hypertension") _tab ("Antivegf") _tab ///
	(`cif_av') _tab (`cif_av_lci') _tab (`cif_av_uci') _tab ///
	(`cif_diff_val') _tab (`cif_diff_lci') _tab (`cif_diff_uci') _tab _n

* Plot
frame cif2 {
	twoway  (rarea F_laser_outcome_lci F_laser_outcome_uci years, color(red%25)) ///
			(rarea F_antivegf_outcome_lci F_antivegf_outcome_uci years, color(blue%25)) ///
			(line F_laser_outcome years, sort lcolor(red)) ///
			(line F_antivegf_outcome years, sort lcolor(blue) lpattern(dash)) ///
			, legend(order(1 "Cataract" 2 "Antivegf") ring(0) cols(1) pos(11) region(lwidth(none))) ///
			title("D", justification(left) size(med)) ///
			yscale(range(0, 35)) ///
			ylabel(0 (5) 35, angle(0) format(%4.1f) labsize(small)) ///
			xlabel(0 (1) 7.3, labsize(small)) ///
			ytitle("Cause-specific CIF (%)", size(medsmall)) ///
			xtitle("Years since index date", size(medsmall)) ///
			graphregion(fcolor(white)) name(nodm_event_hypertension, replace)
}

graph export "$projdir/output/adjcurv_nodm_event_hypertension_cr.tif", as(tif) width(1271) height(762) replace

file close tablech

* Combine eGFR & hypertension cumulative incidence plots 
graph combine dm_egfr_40 dm_event_hypertension ///
              nodm_egfr_40 nodm_event_hypertension, ///
    ycommon xcommon cols(2) ///
    graphregion(color(white)) ///
    name(all_combined, replace)

graph export "$projdir/output/figure_2_revision.tif", ///
    as(tif) width(2542) height(1524) replace
