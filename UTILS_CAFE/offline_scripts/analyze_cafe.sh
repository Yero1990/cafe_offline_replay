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


# Which analysis type are we doing? "data" or "simc"
# User will need to make symlink of the original file (analyze_cafe.sh) to (analyze_cafe_ana_type.sh), where ana_type="simc" or "data"
ana_type=${0##*_}
ana_type=${ana_type%%.sh}  


# Main User Input
runNum=$1     # run number
ana_cut=$2   # CaFe kinematics type, set by user:  "heep_singles", "heep_coin",  "MF", "SRC",
evtNum=$3     # number of events to replay (optional, but will default to all events if none specified)

printDataHelpMsg(){
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./analyze_cafe_${ana_type}.sh <run_number> <ana_cut> <evt_number>"
    echo ""
    echo "<ana_cut> = \"bcm_calib\", \"lumi\", \"optics\", \"heep_singles\", \"heep_coin\", \"MF\" or \"SRC\" "
    echo "<ana_cut> sets analysis cuts based on the type of analysis selected"
    echo ""
    echo "If no <evt_number> specified, defaults to -1 (all events) " 
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

}

printSIMCHelpMsg(){
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Usage:  ./analyze_cafe_${ana_type}.sh <ana_cut> <evt_number>"
    echo ""
    echo "<ana_cut> = \"heep_singles\", \"heep_coin\", \"MF\" or \"SRC\" "
    echo "<ana_cut> sets analysis cuts based on the type of analysis selected" 
    echo ""
    echo "If no <evt_number> specified, defaults to -1 (all events) " 
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 

}


# DATA Help Message / Usage
if [ "${ana_type}" = "data" ]; then
    
    if [ -z "$1" ] || [ -z "$2" ]; then     

	printDataHelpMsg	
	exit 0
	
	# fool-proof make sure only options: bcm_calib, lumi, optics, heep_singles, heep_coin, MF, SRC are used
    elif [ "$ana_cut" == "bcm_calib" ] || [ "$ana_cut" == "lumi" ] || [ "$ana_cut" == "optics" ] || [ "$ana_cut" == "heep_singles" ] ||  [ "$ana_cut" == "heep_coin" ] || [ "$ana_cut" == "MF" ] || [ "$ana_cut" == "SRC" ]; then 
	echo "" 
    else
	printDataHelpMsg	
	exit 0
    fi

    if [ -z "$3" ]; then 
	evtNum=-1
    fi
    
fi

# SIMC Help Message / Usage
if [ "${ana_type}" = "simc" ]; then
    
    if [ -z "$2" ] || [ -z "$3" ]; then     
	
	printSIMCHelpMsg	
	exit 0
	# fool-proof make sure only options: bcm_calib, lumi, optics, heep_singles, heep_coin, MF, SRC are used
    elif [ [ "$ana_cut" == "heep_singles" ] ||  [ "$ana_cut" == "heep_coin" ] || [ "$ana_cut" == "MF" ] || [ "$ana_cut" == "SRC" ] ]; then 
	echo "" 
    else
	printSIMCHelpMsg	
	exit 0
    fi

    if [ -z "$3" ]; then 
	evtNum=-1
    fi
fi


# Default arguments (not required by user as input, unless the user sets as command-line arguments)
daq_mode="coin"
e_arm="SHMS"
hel_flag=0
bcm_type="BCM1"
bcm_thrs=5           # beam current threhsold cut > bcm_thrs [uA]
trig_single="trig2"    # singles trigger type to apply pre-scale factor in FullWeight, i.e. hist->Scale(Ps2_factor) 
trig_coin="trig5"      # coin. trigger type to apply pre-scale factor in FullWeight, i.e., hist->Scale(Ps5_factor)
combine_runs=0

# cafe analysis script
prod_script="UTILS_CAFE/main_analysis.cpp"

# command to run analysis script ( if doing SIMC, all arguments except e_arm, ana_type and ana_cu are irrelevant )
run_cafe="root -l -q -b  \"${prod_script}( ${runNum},    ${evtNum},           
                                    \\\"${daq_mode}\\\",  \\\"${e_arm}\\\",  
                                   \\\"${ana_type}\\\", \\\"${ana_cut}\\\",
                                    ${hel_flag},                        
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_single}\\\", \\\"${trig_coin}\\\", ${combine_runs}
                     )\""



