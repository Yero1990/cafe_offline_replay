#!/bin/bash

#--------------------------------------------
# Author: C. Yero
# Date: Oct 20, 2022
# email: cyero@jlab.org, cyero002@gmail.com
#--------------------------------------------

# Brief: This shell script runs the analyzer on raw replayed
#        .root files (leafs/branches) and outputs .root files
#        containing histogram objects as well as an output .txt
#        file containing a summary of the run and cuts applied. 
#        The user has the option to set/modify analysis cuts,
#        and histogram binning via an external .txt file, as well
#        as applying a scale factor to the data weight for charge
#        normalization, and efficiency corretions, e.g. scale = 1 / (charge * track_eff * live_time *  . . . )
#        By scaling the data by these correction factors, one will be able
#        to make direct comparisons with simulation, as simulation assumes 100% efficiency,
#        and 1 mC of charge.


runNum=$1     # run number
kin_type=$2   # CaFe kinematics type, set by user:  "heep_singles", "heep_coin",  "MF", "SRC", depending on the production type
evtNum=$3     # number of events to replay (optional, but will default to all events if none specified)

# cafe analysis script
prod_script="UTILS_CAFE/main_data_analysis.cpp"

