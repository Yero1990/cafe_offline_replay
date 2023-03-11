Rate Dependence Studies

This directory contains the relevant sub-directories/files used
in the CaFE rate dependence studies. 

The sub-directories definition is:
(the subdirectories phase1-3 below were analyzed with a  >5 uA cut)

phase0: existing CaFe configuration parameters (using imporved pruning tracking algorithm), reference time/ time windows, etc.

phase1: added the T2,T3 TDC raw time window cuts to the param file: PARAM/TRIG/fall22/tcoin_fall22.param 

phase2: added tighter reference time cuts on HMS/SHMS to try and recover more lost coincidences 
(modified: PARAM/SHMS/GEN/fall22/p_reftime_cut_cafe.param  , PARAM/HMS/GEN/fall22/h_reftime_cut_cafe.param  )

phase3: added an event type cut (g.evtyp>=4) to properly select coincidence and EDTM tdc time (for defining the live time)
The following cut was added/modified in the baseAnalyzer.cpp code:  gevtyp>=4, for both data and EDTM TDC time counter
To figure out what event types definition is, see:  DEF-files/CUTS/cafe_cuts.def 
