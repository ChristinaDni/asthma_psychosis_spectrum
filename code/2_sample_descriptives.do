*Sample Descriptives*

use "\dir\revision_03_07_25\dataset\dataset_clean.dta", clear

cd "\dir\revision_03_07_25\observational\results\"

ssc install table1_mc

mark obs
markout obs PE18_24TH_noSF asthma_7 
keep if obs==1


	table1_mc,  by(asthma_7) ///
				vars( /// 
					PE18_24TH_noSF cat %5.1f \ ///
					PE18_24TH_noSF_disORfreq cat %5.1f \ ///
					pliks18_noattSF_current cat %5.1f \ ///
					pliks24TH_noattSF_current cat %5.1f \ ///
					sex cat  %5.1f \ ///
					parity cat %5.1f \ /// 
					mat_education cat %5.1f \ ///   
					fin_problems cat %5.1f \ ///   
					mat_anxiety contn %5.1f \ ///
					mat_depression cat %5.1f \ ///			
					house_smoking cat %5.1f \ ///
					iq_8 contn %5.1f \ ///
					) ///
				nospace onecol missing total(before) test ///
				saving("descriptives_asthma7.xlsx", replace)

			
