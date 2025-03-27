/*==============================================================================
DO FILE NAME:			0_master.do
DATE: 					27/11/2024
AUTHOR:					Ruth Costello 
DESCRIPTION OF FILE:	describes order to run scripts

There are 3 cohorts. Those who have at least one: 
1. anti-vegf injection in HES data
2. photocoagulation in HES data
3. cataracts in HES data
The first injection is during the study period 1st January 2013 - 29th February 2020

==============================================================================*/

do "$dodir/1_define_study_population.do"

do "$dodir/2_define_variables.do"

do "$dodir/2_define_outcomes.do"

do "$dodir/3_descriptives.do"

do "$dodir/3_main_analysis.do"

do "$dodir/3_sensitivity_analysis.do"

do "$dodir/3_stratified_analysis.do"

do "$dodir/3_forest_plots.do"
