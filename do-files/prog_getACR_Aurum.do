/*=========================================================================
DO FILE NAME:	prog_getACR.do

AUTHOR:					Paris Baptiste

VERSION:				v1
DATE VERSION CREATED:	21-Sept-2020					
					
DATABASE:				CPRD July 2019 build

DESCRIPTION OF FILE: 
	Extracts ACR test results.
	Flags if proteinuria
	
	Arguements (options) required:
		* testfile			// path and name of file containing test result 
		* testfilesnum 		// number of test files to loop through
		* savefile			// string containing name of file to save

									
HOW TO USE: e.g.

run prog_getACR.do

prog_getACR, ///
	testfile("$RawData\Test_extract_CVD") ///
	testfilesnum(5) ///
	savefile("$DataDerived\cases-Proteinuria") ///

*=========================================================================*/


/*******************************************************************************
#>> Define program
*******************************************************************************/
capture program drop prog_getACR_Aurum
program define prog_getACR_Aurum

syntax, obsfile(string) obsfilenum(integer) proteinuria_codelist(string) measure(string) savefile(string) 
	
* testfile			// path and name of file containing test result extract files (exclude the underscore and the number of the file)
* testfilesnum 		// number of test files to loop through
* codelist			// list of medcodes that are likely to be used for ACR test results
* savefile			// string containing name of file to save

noi di
noi di in yellow _dup(15) "*"
noi di in yellow "Identify ACR test results and identify proteinuria"
noi di in yellow _dup(15) "*"


qui{
/*******************************************************************************
================================================================================
#A. EXTRACT AND CLEAN ACR RESULTS
================================================================================
*******************************************************************************/

	/*******************************************************************************
	#A1. Identify test records for ACR results.
	*******************************************************************************/
	display in red "*******************test file number: 1*******************"

	use "`obsfile'_1", clear
	
	merge m:1 ptid using "Z:\GPRD_GOLD\Paris\ACEi_ARB\Aurum\Extract_from_start_2001\CPRD extract HES linked\eligible periods\Drug_issue_CVD_pts.dta", keepusing(ptid) keep(match) nogen
	
	merge m:1 medcodeid using "`proteinuria_codelist'", nogen keep(match) force
	save `savefile', replace
	
	

	/*******************************************************************************
	#A2. Loop through subsequent (from 2 onwards) separate test extract files in 
		turn and append the results to the first extract file saved in #1
	*******************************************************************************/
	forvalues n=2/`obsfilenum' {
		display in red "*******************test file number: `n'*******************"

		use "`obsfile'_`n'", clear
		merge m:1 ptid using "Z:\GPRD_GOLD\Paris\ACEi_ARB\Aurum\Extract_from_start_2001\CPRD extract HES linked\eligible periods\Drug_issue_CVD_pts.dta", keepusing(ptid) keep(match) nogen
		* keep only relevant medcodes (urine albumin:creatinine ratio, albumin/creatinine 	ratio, urine microalbumin:creatinine ratio)
		merge m:1 medcodeid using "`proteinuria_codelist'", nogen keep(match) force

		* add the file containing records for the specified comorbidity
		* to make one file containing all specified comorbidiy records for the
		* clinical extract specified
		append using "`savefile'"
		
		* save
		save "`savefile'", replace
	}
	
	
	
	
	/*******************************************************************************
	#A3. Drop unnecessary vars and label variables.
	*******************************************************************************/	
	* drop unuseful vars
	*drop constype consid staffid data8
	 
	*rename variables and add labels
	*rename data1 operator		// e.g. 1"<" 2"<="
	capture destring value, replace
	rename value `measure' 			// ACR result
	rename numunitid unit 			// unit of measure
	*rename data4 qualifier		// e.g. 7"very high" 
	rename numrangelow rangeFrom 		//"normal range from"
	rename numrangehigh rangeTo		//"normal range to"
	*rename data7 rangeBasis		//"normal range basis"
	
	*label variable operator "operator: eg 1 is <, 2 is <="	
	label variable `measure' "`measure': `measure' result"
	label variable unit "unit of measure"	
	*label variable qualifier "qualifier e.g. 7-very high" 
	label variable rangeFrom "rangeFrom: normal range from"
	label variable rangeTo "rangeTo: normal range to"
	*label variable rangeBasis "rangeBasis: normal range basis"
	
	*duplicates drop 
	
	*Label values (info from enttype 469 look-up file
	*#delimit;
		*label define UNIT 0 "[0] missing" 1 "[1] %" 15 "[15] /day" 21 "[21] /mL" 26 "[26] 1" 41 "[41] d" 47 "[47] fL" 49 "[49] g" 50 "[50] g(creat)" 55 "[55] g/d" 57 "[57] g/L" 58 "[58] h" 60 "[60] IU/d" 74 "[74] m" 80 "[80] mg" 81 "[81] mg/d" 82 "[82] mg/dL" 83 "[83] mg/L" 84 "[84] mg/m3" 85 "[85] mg/min" 86 "[86] mg/mmol" 87 "[87] min" 90 "[90] mL/min" 92 "[92] mm" 94 "[94] mm/h" 96 "[96] mmol/L" 97 "[97] mmol/mol" 99 "[99] mol/L"  104 "[104] ng" 105 "[105] ng/L" 106 "[106] ng/mL" 108 "[108] nm" 110 "[110] nmol/L"  122 "[122] rad" 124 "[124] s" 126 "[126] u" 133 "[133] ug/L"  138 "[138] umol" 140 "[140] umol/g(creat)" 142 "[142] umol/L" 147 "[147] week" 149 "[149] MicroU/L" 151 "[151] ratio" 154 "[154] L" 156 "[156] mmol" 160 "[160] mu/L" 161 "[161] 1/1" 165 "[165] IU/dL" 169 "[169] mg/mmol(creat)" 171 "[171] mmol/mmol(creat)" 178 "[178] u/mL" 185 "[185] L/L" 187 "[187] mmol/mmol" 188 "[188] nmol/g" 191 "[191] nmol/mmol" 193 "[193] uU/spec" 194 "[194] g/mol" 195 "[195] titre" 196 "[196] mg/12hrs" 210 "[210] #/HPF" 216 "[216] /g(Hb)" 222 "[222] fg/L" 225 "[225] g/kg" 239 "[239] Mcell/L" 241 "[241] mg/24h" 242 "[242] mg/g(dry wt)" 243 "[243] mg/g(wet wt)" 247 "[247] mg/mg" 254 "[254] mmol/g(wet wt)" 262 "[262] nmol/mmol(creat)" 263 "[263] nmol/mol(creat)" 273 "[273] ug/g(wet wt)" 275 "[275] ug/mmol(creat)" 282 "[282] umol/h/L" 296 "[296] mg/g" 299 "[299] N/A" 316 "[316] U/mmol" 317 "[317] ug/mmol" 322 "[322] mg/mL" 330 "[330] Ratio";
		*label values unit UNIT;

		*label define OPR 0 "[0] missing" 1 "[1] <" 2 "[2] <=" 3 "[3] =" 4 "[4] >" 5 "[5] >=" 6 "[6] ~";
		*label values operator OPR;

		*label define QUALIFIER 0 "[0] missing" 7 "[7] Very High" 8 "[8] High" 9 "[9] Normal" 10 "[10] Low" 11 "[11] Very Low" 12 "[12] Abnormal" 13 "[13] Potential Abnormal" 14 "[14] Outside ref range" 35 "[35] Above high reference limit" 38 "[38] Potentially abnormal" 39 "[39] Normal" 40 "[40] Low" 41 "[41] High" 44 "[44] Abnormal";
		*label values qualifier QUALIFIER;

		*label define POP 0 "[0] missing" 1 "[1] Age based" 3 "[3] Age and sex based" 5 "[5] Generic normal range" 6 "[6] Other (e.g.specific disease based)" 9 "[9] Sex based";
		*label values rangeBasis POP;
	*#delimit cr;
	
	
	
	/*******************************************************************************
	#A4. Drop any duplicate records
		Drop records with missing dates or ACR results
	*******************************************************************************/	
	duplicates drop

	* drop if eventdate missing 
	* but check if sysdate available and replace missing eventdate with sysdate if available
	replace obsdate=enterdate if (obsdate=="" & enterdate!="")
	gen obsdate1=date(obsdate,"DMY")
	format obsdate1 %td
	
	drop obsdate 
	rename obsdate1 obsdate
	drop if obsdate==.
	
	
	* drop if ACR value is missing or zero
	drop if `measure'==0 
	drop if `measure'==.

	
	
	/*******************************************************************************
	#A6. Drop records with ACR values that are very low or very high
	*******************************************************************************/
	if "`measure'"=="ACR" {
		* drop improbable values for ACR i.e. <0 or >3000
		gen improbable=0
		recode improbable 0=1 if ACR<0 | ACR>2000 | ACR==0

		drop if improbable==1
		drop improbable	
	}
	
	if "`measure'"=="potassium" {
		* drop improbable values for ACR i.e. <0 or >3000
		gen improbable=0
		recode improbable 0=1 if potassium<0 | potassium>10 | potassium==0

		drop if improbable==1
		drop improbable	
	}
	
	/*******************************************************************************
	#A7. Add notes and labels to the data
	*******************************************************************************/
	notes: prog_getProteinuria.do / TS
	save "`savefile'", replace

	
	
	
/*******************************************************************************
================================================================================
#B. Flag is ACR>3
================================================================================
*******************************************************************************/
if "`measure'"=="ACR" {
	gen proteinuria=1 if ACR>3
	replace proteinuria=0 if proteinuria==.
	display in red "*******************final dataset*******************"	

	* save	
	label data "proteinuria records from CPRD"
	notes: prog_getProteinuria.do / TS
	save "`savefile'", replace
}
}/*end of quietly*/

end


