#! /bin/bash


#--------------------------------------------
# Author: C. Yero
# Date: Oct 20, 2022
# email: cyero@jlab.org, cyero002@gmail.com
#--------------------------------------------

# Display help output if no argumenr specified
if [  $# -eq 0 ]; then
    echo "" 
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:"
    echo ""
    echo "Brief: This shell script replays the raw (.dat) files "
    echo "       and outputs the raw (.root) files. The leaf "
    echo "       variables to be written are specified under "
    echo "       DEF-files/cafe_prod.def and may be modified. "
    echo ""
    echo "------------------------------------"
    echo "Usage 1):  ./replay_cafe.sh run evt "
    echo "------------------------------------"
    echo ""
    echo "run: run number"
    echo ""
    echo "evt: event number defaults to -1 (all events), "
    echo "unless explicitly specified as 3rd argument"
    echo ""
    echo "example 1: ./replay_cafe.sh 3288 100000"
    echo ""
    echo "-------------------------------------------"
    echo "Usage 2):  ./replay_cafe.sh target kin evt "
    echo "-------------------------------------------"
    echo ""
    echo "target:  dummy, h, d2, be9, b10, b11, c12, ca40, ca48, fe54 "
    echo ""
    echo "kin: singles, heep_coin, MF or SRC "
    echo ""
    echo "evt: event number defaults to -1 (all events), "
    echo "unless explicitly specified as 3rd argument"
    echo ""
    echo "example 2: ./replay_cafe.sh ca48 MF 350000"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:" 
    echo "" 
    exit 0  
fi


run=$1

# replay script
replay_script="SCRIPTS/COIN/PRODUCTION/replay_cafe.C" 

# check if 1st argument is an integer (i.e., run number, else attempt to read from runlist)
if [ $1 -eq $1 ]; then

    evt=$2

    # check if event number is specified
    if [ -z $evt ]; then
	evt=-1
	echo "No event number spedified, defaulting to evt=${evt} (all events)"
    fi

    # hcana command
    run_hcana="./hcana -q \"${replay_script}(${run}, ${evt}, \\\"prod\\\")\""

    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    echo "" 
    date
    echo ""
    echo ""
    echo "Running HCANA CaFe Replay on the run ${run}:"
    echo " -> SCRIPT:  ${replay_script}"
    echo " -> RUN:     ${run}"
    echo " -> NEVENTS: ${evt}"
    echo " -> COMMAND: ${run_hcana}"
    echo ""
    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
    
    sleep 2
    eval ${run_hcana}

    
else

    # order of arguments now is as follows: target, kin, evt (run number excluded)
    # read 1st and 2nd arguments
    target=$1
    kin=$2
    
    evt=$3
    # check if event number is specified
    if [ -z $evt ]; then
	evt=-1
	echo "No event number spedified, defaulting to evt=${evt} (all events)"
    fi
    
    # runlist
    filename="UTILS_CAFE/runlist/${target}_${kin}.txt"

    for run in $(cat $filename) ; do    

	# hcana command
	run_hcana="./hcana -q \"${replay_script}(${run}, ${evt}, \\\"prod\\\")\""
    
	{
	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    echo "" 
	    date
	    echo ""
	    echo ""
	    echo "Running HCANA CaFe Replay on the run ${runNum}:"
	    echo " -> RUNLIST: ${filename}"
	    echo " -> SCRIPT:  ${replay_script}"
	    echo " -> RUN:     ${run}"
	    echo " -> NEVENTS: ${evt}"
	    echo " -> COMMAND: ${run_hcana}"
	    echo ""
	    echo ":=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:="
	    
	    sleep 2
	    eval ${run_hcana}
	}


    done
fi







  
