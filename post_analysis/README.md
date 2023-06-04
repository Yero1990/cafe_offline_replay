# CaFe Post-Analysis

## Introduction
This directory stores the relevant CaFe analysis
output files (.root, .csv) and post-analysis utility scripts for making production plots. 

### Directory Structure:
* ` analyzed_files/individual/passN/`<br> 

	This directory contains the analyzed individual run-by-run ROOTfiles and REPORT_FILES of a particular (`target`, `kinematics`), with a generic filename: <br>

        cafe_prod_<target>_<kin>_<run>_-1_histos.root 
        cafe_prod_<target>_<kin>_report_<run>_-1.txt        


* ` analyzed_files/combined/passN/`<br> 

	This directory contains the analyzed ROOTfiles summed over all runs of a particular (`target`, `kinematics`), with a generic filename: <br>

        `cafe_prod_<target>_<kin>_combined.root` 

	Each individual (or combined) ROOTfile consist of sub-directory structure with different levels of histogram categories as follows (and may be subject to new additions):
	* `quality_plots/`: contains 1) different levels of cuts to study their effects on the yield, and 2) fits to calibration check histograms for quality monitoring 
	
	* `pid_plots/`: particle identification (pid)-related histograms with all analysis cuts + main coincidence time cut
	
	* `kin_plots/`: kinematic-related histograms with all analysis cuts + main coincidence time cut
	
	* `accp_plots/`: detector acceptance-related histograms with all analysis cuts + main coincidence time cut
	
	* `rand_plots/`: selected histograms with all analysis cuts + random coincidence time selection (for accidental background subtraction)
	
	* `randSub_plots/`: selected histograms with randoms-subtracted analysis cuts (these histograms represent the true (or real) background subtracted yield that will be used in post-analysis 
	
	**NOTE:** 
	
	1) Calcium (Ca-48) target needs to be corrected for (H, C) contamination on a run-by-run basis (< ~3%) + Ca-40 impurity (~10 %), 
	
	2) Boron Carbide isotopes, (B4C)-10, (B4C)-11 targets need to be Carbon (C12)-subtracted

 These corrections are done separately in the script: `scripts/make_final_summary.py`
	
* 	`summary_files/passN/` <br>
	This directory contains report summary files over all runs of a particular (`target`, `kinematics`), with a generic filename: <br>

 `cafe_prod_<target>_<kin>_summary.csv` <br>

	Each .csv file consists of multi-column data on a run-by-run basis, and which ultimately need to be combined for taking ratios. For example, the **charge** column needs to be added to get total charge, the **real_Yield** column needs to be added to get total real counts (either from MF or SRC), and the efficiency columns need to be averaged over all runs assuming they are ~stable on a run-basis
	
* 	`scripts/` <br>
This directory contains relevant CaFe-specific utility scripts for histogram comparisons and summary files computations.

* 	`special_studies/` <br>
This directory contains sub-directories with specific special studies done for the
CaFe post-replay analysis stage. These include target impurity subtraction, Ca48 contamination corrections, rate-dependence studies, H(e,e'p) checks and alignments, systematics studies, etc.
