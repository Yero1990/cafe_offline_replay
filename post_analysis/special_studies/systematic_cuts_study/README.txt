-----------------------
systematics_cuts_study
-----------------------

This directory contains the relevant information
concerning the CaFe systematics studies done by
varying a given set of kinematical cuts about a
randomly-discributed gaussian distribution with
mean = central cut value and a spread of 2sigma

------------------------------
Relevant Directory Structure
------------------------------

scripts/:   relevant systematic analysis scripts

1) genRandCuts.py     -> generates numerical file with N random
                       	 gaussian-distributed combination of cuts
			 as well as the plots with the histogrammed
			 numerical values
			 
2) fit_systematics.py -> fits the numerical data output (integrated yield)
                         results from varying the cuts N times. The
			 integrated yield counts represent the counts
			 obtained for every combination of analysis cuts
			 The fits are done for: invidivual (target,kin),
			 single and double ratios

input/: contains the input numerical cuts file to be read in by
        the analysis code. The identical analysis is performed for
	every entry in this file as part of the systematic cuts study


output/: contains the output numerical files with the integrated yield results
         where each entry represents the results obtained from applying the
	 corresponding combination of cuts from the input/ file
	 
         
