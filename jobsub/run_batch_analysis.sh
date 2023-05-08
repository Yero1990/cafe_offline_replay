#! /bin/bash

### Stephen Kay, University of Regina
### 03/03/21
### stephen.kay@uregina.ca
### A batch submission script based on an earlier version by Richard Trotta, Catholic University of America
### SJDK - 06/01/22 - Updated to use swif2 system, also cleaned up the script and added some more comments

# usage: ./run_batch_template.sh ca40_MF.txt

echo "Running as ${USER}"
#RunList=$1

# user input
targ=$1
ana_cut=$2
RunList="${targ}_${ana_cut}.txt"

echo "RunList: ${RunList}"
inputFile="/w/hallc-scshelf2102/c-cafe-2022/cyero/cafe_offline_replay/UTILS_CAFE/runlist/"

if [[ -z "$1" ]]; then
    echo ""
    echo "-----------------------------------------------------"
    echo ""
    echo "Usage: "
    echo "./run_batch_analysis.sh <targ> <ana_cut> MAXEVENTS"
    echo ""
    echo "<targ> : target type (LH2, LF2, Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197) "
    echo "<ana_cut> : analysis type cut (heep_coin, heep_singles, MF, SRC, optics)"
    echo ""
    echo "A runlist will be read, based on the input, from the directory located at: "
    echo "${inputFile}<targ>_<ana_cut>.txt"
    echo ""
    echo "MAXEVENTS : maximum number of events to be replayed"
    echo "defaults to all events (-1) if no argumnet provided "
    echo ""
    echo "-----------------------------------------------------" 
    echo ""
    exit 2
fi
# Check if an argument was provided, if not assume -1, if yes, this is max events
if [[ $3 -eq "" ]]; then
    MAXEVENTS=-1
else
    MAXEVENTS=$3
fi

# 15/02/22 - SJDK - Added the swif2 workflow as a variable you can specify here
#Workflow="cafe_aug08_${USER}" # Change this as desired
Workflow="cafe_analysis_${USER}" # Change this as desired

# Input run numbers, this just points to a file which is a list of run numbers, one number per line
inputFile="${inputFile}${RunList}"

if test -f "$inputFile"; then
    echo "Reading the input file: ${inputFile}" 
else
    echo "${inputFile} does NOT EXIST !"
    exit 0
fi
# Output batch job text file, this is the script that is submitted as part of the job, 
# change the name of this as you want. Preferably, you should have different script name for each job so that you don't get any overwriting weirdness##
#batch="${USER}_${runNum}_Job.txt"
batch="${USER}_${targ}_${ana_cut}_Job.txt"

cp /dev/null ${batch}
# Creation of batch script for submission, the script is just a series of commands
# We add these to the file (${batch}) via a series of piped echo commands
echo "PROJECT: c-comm2017" >> ${batch} # Or whatever your project is!
echo "TRACK: analysis" >> ${batch} ## Use this track for production running
#echo "TRACK: debug" >> ${batch} ### Use this track for testing, higher priority
echo "JOBNAME: cafe_${targ}_${ana_cut}" >> ${batch} ## Change to be more specific if you want
# Request double the tape file size in space, for trunctuated replays edit down as needed
# Note, unless this is set typically replays will produce broken root files
echo "DISK_SPACE: 5 GB" >> ${batch}
echo "MEMORY: 10000 MB" >> ${batch} 
echo "CPU: 1" >> ${batch} ### hcana is single core, setting CPU higher will lower priority and gain you nothing!
echo "COMMAND:/w/hallc-scshelf2102/c-cafe-2022/cyero/cafe_offline_replay/analyze_cafe_data.sh ${targ} ${ana_cut} ${MAXEVENTS}"  >> ${batch}
# simulation script (will need to make alternate submittion script for simulation)
#echo "COMMAND:/w/hallc-scshelf2102/c-cafe-2022/cyero/hallc_simulations/simulate.py" >> ${batch}
echo "MAIL: ${USER}@jlab.org" >> ${batch}
echo "Submitting batch"
# swif2 is now used for job submission, we use our old jsub style scripts. The argument set to "LTSep" currently is the workflow. Change this if you want.
eval "swif2 add-jsub ${Workflow} -script ${batch} 2>/dev/null" # Swif2 job submission, uses old jsub scripts
echo " "
sleep 2
# Delete the script we just submitted as a batch job, this stops this folder getting clogged
rm ${batch}



