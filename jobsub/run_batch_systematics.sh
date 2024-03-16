#! /bin/bash

# run analysis code to do systematics study on cafe
# usage: ./run_batch_systematics.sh target kin

echo "Running as ${USER}"

# user input
targ=$1
ana_cut=$2

if [[ -z "$1" ]]; then
    echo ""
    echo "-----------------------------------------------------"
    echo ""
    echo "Usage: "
    echo "./run_batch_systematics.sh <targ> <ana_cut> "
    echo ""
    echo "<targ> : target type (LD2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197) "
    echo "<ana_cut> : analysis type cut (MF, SRC)"
    echo ""
    echo "-----------------------------------------------------" 
    echo ""
    exit 2
fi


Workflow="cafe_systematics_${USER}" # Change this as desired


# Output batch job text file, this is the script that is submitted as part of the job, 
batch="${USER}_${targ}_${ana_cut}_Job.txt"

cp /dev/null ${batch}
# Creation of batch script for submission, the script is just a series of commands
# We add these to the file (${batch}) via a series of piped echo commands
echo "PROJECT: c-comm2017" >> ${batch} # Or whatever your project is!
echo "TRACK: analysis" >> ${batch} ## Use this track for production running
#echo "TRACK: debug" >> ${batch} ### Use this track for testing, higher priority
echo "JOBNAME: cafe_${targ}_${ana_cut}_syst" >> ${batch} ## Change to be more specific if you want
# Request double the tape file size in space, for trunctuated replays edit down as needed
# Note, unless this is set typically replays will produce broken root files
echo "DISK_SPACE: 5 GB" >> ${batch}
echo "TIME: 100:00:00" >> ${batch}
echo "MEMORY: 10000 MB" >> ${batch} 
echo "CPU: 8" >> ${batch} ### hcana is single core, setting CPU higher will lower priority and gain you nothing!
echo "COMMAND:/w/hallc-scshelf2102/c-cafe-2022/cyero/cafe_offline_replay/analyze_cafe_systematics.sh ${targ} ${ana_cut}"  >> ${batch}
# simulation script (will need to make alternate submittion script for simulation)
#echo "COMMAND:/w/hallc-scshelf2102/c-cafe-2022/cyero/hallc_simulations/simulate.py" >> ${batch}
echo "MAIL: ${USER}@jlab.org" >> ${batch}
echo "Submitting batch"
# swif2 is now used for job submission, we use our old jsub style scripts. The argument set to "LTSep" currently is the workflow. Change this if you want.
eval "swif2 add-jsub ${Workflow} -script ${batch} 2>/dev/null" # Swif2 job submission, uses old jsub scripts
echo "-----> ${batch}"
echo " "
sleep 2
# Delete the script we just submitted as a batch job, this stops this folder getting clogged
rm ${batch}



