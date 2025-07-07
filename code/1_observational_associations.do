
* Association analyses

use "\dir\revision_03_07_25\dataset\dataset_clean.dta", clear

cd "\dir\revision_03_07_25\observational\results"

*Crude analyses between all exposures-outcomes

global obs_expvars asthma_7
global outvars PE18_24TH_noSF PE18_24TH_noSF_disORfreq Ppliks18_noattSF_current pliks24TH_noattSF_current 
global obs_covars sex parity mat_education fin_problems mat_anxiety mat_depression house_smoking iq_8

*Unadjusted analyses

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI  ///
using "\dir\revision_03_07_25\observational\results\crude_exposure_outcome.dta", replace

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in PE18_24TH_noSF PE18_24TH_noSF_disORfreq pliks18_noattSF_current pliks24TH_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' 
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}


postclose `memhold'

*PEs since the age of 12
*Exposure-outcome-covariate data 
*Unadjusted 

mark obs
markout obs PE18_24TH_noSF PE18_24TH_noSF_disORfreq asthma_7 sex parity mat_education fin_problems mat_anxiety mat_depression house_smoking iq_8
keep if obs==1

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI  ///
using "\dir\revision_03_07_25\observational\results\psychotic_experiences_since12_associations.dta", replace

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in PE18_24TH_noSF PE18_24TH_noSF_disORfreq {
		local i = `i' + 1		
		logistic `outcome' `exposure' 
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Unadjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}

*Exposure-outcome-covariate data
*Adjusted

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in PE18_24TH_noSF PE18_24TH_noSF_disORfreq {
		local i = `i' + 1		
		logistic `outcome' `exposure' i.sex i.parity i.mat_education i.fin_problems i.mat_depression mat_anxiety iq_8 i.house_smoking
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Adjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}


postclose `memhold'

*PEs past six months at 18
*Exposure-outcome-covariate data 
use "\dir\revision_03_07_25\dataset\dataset_clean.dta", clear

*Unadjusted 

mark obs_18
markout obs_18 pliks18_noattSF_current asthma_7 sex parity mat_education fin_problems mat_anxiety mat_depression house_smoking iq_8
keep if obs_18==1

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI  ///
using "\dir\revision_03_07_25\observational\results\psychotic_experiences_past6months_18_associations.dta", replace

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in pliks18_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' 
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Unadjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}

*Exposure-outcome-covariate data
*Adjusted

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in pliks18_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' i.sex i.parity i.mat_education i.fin_problems i.mat_depression mat_anxiety iq_8 i.house_smoking
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Adjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}


postclose `memhold'


*PEs past six months at 24
*Exposure-outcome-covariate data 
use "\dir\revision_03_07_25\dataset\dataset_clean.dta", clear

*Unadjusted 

mark obs_24
markout obs_24 pliks24TH_noattSF_current asthma_7 sex parity mat_education fin_problems mat_anxiety mat_depression house_smoking iq_8
keep if obs_24==1

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI  ///
using "\dir\revision_03_07_25\observational\results\psychotic_experiences_past6months_24_associations.dta", replace

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in pliks24TH_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' 
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Unadjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}

*Exposure-outcome-covariate data
*Adjusted

foreach exposure in asthma_7 {
	local i=0
	foreach outcome in pliks24TH_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' i.sex i.parity i.mat_education i.fin_problems i.mat_depression mat_anxiety iq_8 i.house_smoking
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Complete data on exposure outcome covariates- Adjusted") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}


postclose `memhold'



