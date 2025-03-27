/******************************************************************************* 

NAME: 			make_episodes.do

PURPOSE: 		defines a program, which can be called in other dofiles 
				the program makes treatment episodes for a single prescription type 
				user can specify how overlaps should be handled: 
				
					1) ignored 
					2) added
					3) added, but total number of days added capped 
				
				flexible input for allowable gap 
				
AUTHOR:			anna schultze, based on code from Angel Wong
				original code logic from statalist, here: https://www.statalist.org/forums/forum/general-stata-discussion/general/1349596-overlapping-time-intervals-how-to-combine
				
DATA USED: 		input dataset should only contain drug you want to make episodes for
				each dataset should contain only the variables used in the program call 

EDITS: 			n/a

*******************************************************************************/ 

********************************************************************************

/* Explanatory Notes 
defines program make_episodes
syntax is : 
	patid - your patid variable 
	start_date - your start date variable 
	end_date - your end date variable 
	allowable_gap - the allowable gap between an end and start date  
	overlap - how to handle overlaps, option "ignore", "add", "cap"
	overlap_cap - OPTIONAL, if overlap == cap, what's the maximum added duration? 
				  the number of maximum "carried over" days from one prescription to the next in days  
*/ 

********************************************************************************
cap prog drop make_episodes 
program define make_episodes
syntax, patid(varname) start_date(varname) end_date(varname) allowable_gap(string) overlap(string) [overlap_cap(string)]

	di "NOTE: initiating make_episodes"

	* 1) CREATE VARIABLES NEEDED================================================
				
		gen time_start = `start_date'
		gen time_end = `end_date'
		
		drop `start_date' `end_date'
		
		bysort `patid' (time_start): gen episode = _n 
		
		preserve 
		duplicates drop patid, force
		count 
		di "r(N) individuals in this dataset after processing"
		restore 
		
	* 2) HANDLE OVERLAPS BETWEEN PRESCRIPTIONS==================================	
		* options for how overlapping prescriptions should be handled 
		
		* add all 'spare' days at a new prescription to the end of the new one 
		if "`overlap'" == "add" {
			
			di "NOTE: overlap option ADD - prescription overlaps are ADDED"
			
			* calculate the gap between a new prescription and when the previous one ended
			bysort `patid': gen rxdiff = (time_end[_n - 1] - time_start) 
			
			* need to generate a cumulative duration variable 
			* this should sum the extra number of days resulting from overlap, so we only want to sum when the difference is +
			* however, we sometimes want to sum when the difference is - if the resulting sum is still positive 
			* this is because the overlap is "carried over" with each prescription 
			* therefore create a cumulative sum, which sums the prescription gap ONLY if that sum is greater than zero
			bysort `patid': gen cum_overlap = 0
			bysort `patid': replace cum_overlap = max((cum_overlap[_n-1] + rxdiff),0) 
			
			* print the distribtion for reference 
			di "CHECK: Distribution of overlapping days that are added to prescriptions"
			summ cum_overlap, d 
			
			* calculate new end dates
			replace time_end = time_end + cum_overlap
			format time_start time_end %td 
			drop cum_overlap rxdiff 
			
		}
		
		* as above, but cap the days carried over to the value specified in cap 
		else if "`overlap'" == "cap" {
			
			capture assert (`overlap_cap' != . & `overlap_cap' > 0) 
			
			if _rc == 9 {
				    
					di "cap on overlap days (overlap_cap) not specified in call, cannot proceed"
				    
			} 
			
			else if _rc == 0 {
			    
			di "NOTE: overlap option CAP - adding prescription overlaps, MAX `overlap_cap' days"
				
			* calculate the gap between a new prescription and when the previous one ended
			bysort `patid': gen rxdiff = (time_end[_n - 1] - time_start) 
			
			* as above, but cap the cumulative overlap at + 90
			bysort `patid': gen cum_overlap2 = 0
			bysort `patid': replace cum_overlap2 = max((min(cum_overlap2[_n-1] + rxdiff,`overlap_cap')),0) if rxdiff !=.
			
			* print the distribtion for reference 
			di "CHECK: Distribution of overlapping days that are added to prescriptions"
			summ cum_overlap2, d 
			
			* calculate new end dates
			replace time_end = time_end + cum_overlap2
			format time_start time_end %td 
			drop cum_overlap2 rxdiff 
			
			}
		
		}
		
		* ignore, prescription assumed to end at issue of a new one 
		else if "`overlap'" == "ignore" {
			
			di "NOTE: overlap option IGNORE - ignoring overlaps between prescriptions"
			
		}
		
	* 3) DATA MANAGEMENT========================================================		
		* allowable gap needs to be added and data needs to be long 
		
		* extend the duration with the allowable gap 
		gen gap = real("`allowable_gap'")
		replace time_end = time_end + gap
		
		* required format has a row for each 'change', ie start and end date
		reshape long time, i(`patid' episode) j(event_type) string
		format time %td 
		
		* generate a variable to ensure that if a start and end date overlap exactly, the start date comes before the end date 
		* this ensures that these kinds of prescriptions are counted as a single episode rather than two separate ones 
		gen event_order = 1 if event_type == "_start"
		replace event_order = 2 if event_type == "_end"

		sort `patid' time event_order
		
	* 4) MAKE TREATMENT EPISODES================================================
	
		* for each person: 
		* count - cumulatively - the number of starts, and subtract the cumulative number of end dates
		* for the first start date this will be one. If there is no overlap then this will be followed be an end-date, and count 0 
		* when an exposure period 'ends' this count will always be represented by a zero 
		* do not sort on episode. If you do this, start and end dates which are the same will contribute to different episodes 
		bysort `patid' (time): gen cont_exposure = sum(event_type == "_start") - sum(event_type == "_end")
		
		* all counts of 1 or greater represent continous exposure 
		replace cont_exposure = 1 if cont_exposure > 1
		
		* count the number of distinct 'episodes' with PID
		by `patid' (time): gen treatment_episode = 1 if cont_exposure == 1 & cont_exposure[_n-1] != 1
		by `patid' (time): replace treatment_episode = sum(treatment_episode)
		
		sort `patid' treatment_episode time 
		
		* keep the first and last row within each episode (to reshape back to wide with a start and end date)
		by `patid' treatment_episode (time): keep if _n == 1 | _n == _N
		 
		drop episode cont_exposure event_order 
		reshape wide time, i(`patid' treatment_episode) j(event_type) string
		
	* 4) FORMAT DATA FOR SAVING=================================================
		rename time_* * 	
		order start, before(end)
		
		gen actual_end = end - gap 
		format actual_end %td 
		
		drop gap 
		
		preserve 
		duplicates drop patid, force
		count 
		di "r(N) individuals in this dataset after processing"
		restore 
	
end
