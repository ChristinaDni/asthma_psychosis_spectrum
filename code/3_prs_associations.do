*PRS analyses

use "\dir\revision_03_07_25\dataset\dataset_clean.dta", clear

cd "\dir\revision_03_07_25\prs\results\"

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI AUC ///
using "\dir\revision_03_07_25\prs\results\asthma_validation.dta", replace

*This bit of script updated to add AUC estimate for each model

forval k=1/13 {
foreach exposure in zscore_asthma_prs_child_S`k' {
	local i=0
	foreach outcome in asthma_7 {
		local i = `i' + 1		
		logistic `outcome' `exposure' pc1 pc2 pc3 pc4 pc5 ///
		pc6 pc7 pc8 pc9 pc10 sex
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		lroc, nograph
		local auc = r(area) 
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Adjusted for PCs and sex") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")") ///
		(`auc')
	}
}
}

postclose `memhold'

capture postutil close
tempname memhold

postfile `memhold' _outord _modelord str50 exposure str150 outcome str50 model modn modn_wout or lci uci p str30 OR_CI  ///
using "\dir\revision_03_07_25\prs\results\asthma_association.dta", replace

forval k=1/13 {
foreach exposure in zscore_asthma_prs_child_S`k' {
	local i=0
	foreach outcome in PE18_24TH_noSF PE18_24TH_noSF_disORfreq pliks18_noattSF_current pliks24TH_noattSF_current {
		local i = `i' + 1		
		logistic `outcome' `exposure' pc1 pc2 pc3 pc4 pc5 ///
		pc6 pc7 pc8 pc9 pc10 sex
		local or  = r(table)[1,1]
		local lci = r(table)[5,1]
		local uci = r(table)[6,1]
		local p= r(table)[4,1]
		local modn = e(N)
		count if e(sample) & `outcome' == 1
		local modn_wout = r(N)
		post `memhold' (`i') (1) ("`exposure'") ("`outcome'") ("Adjusted for PCs and sex") (`modn') (`modn_wout') ///
		(`or') (`lci') (`uci') (`p') ///
		(strofreal(`or', "%5.2f") + " (" + strofreal(`lci',"%5.2f") + "-" + strofreal(`uci',"%5.2f") + ")")             
	}
}
}

postclose `memhold'

