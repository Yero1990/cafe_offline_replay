#! /bin/bash

### Stephen Kay, University of Regina
### 03/03/21
### stephen.kay@uregina.ca
### A batch submission script based on an earlier version by Richard Trotta, Catholic University of America
### SJDK - 06/01/22 - Updated to use swif2 system, also cleaned up the script and added some more comments

# usage: ./run_batch_template.sh ca40_MF.txt

echo "Running as ${USER}"
RunList=$1

inputFile="/w/hallc-scshelf2102/c-cafe-2022/cyero/cafe_offline_replay/UTILS_CAFE/runlist/"

if [[ -z "$1" ]]; then
    echo ""
    echo "-----------------------------------------------------"
    echo ""
    echo "Usage: "
    echo "./run_batch_template.sh runlist.txt MAXEVENTS"
    echo ""
    echo "runlist.txt : a simple text file with a list of "
    echo "run numbers (one run per line), assumed to be on:"
    echo "${inputFile}"
    echo ""
    echo "MAXEVENTS : maximum number of events to be replayed"
    echo "defaults to all events (-1) if no argumnet provided "
    echo ""
    echo "-----------------------------------------------------" 
    echo ""
    exit 2
fi
# Check if an argument was provided, if not assume -1, if yes, this is max events
if [[ $2 -eq "" ]]; then
    MAXEVENTS=-1
else
    MAXEVENTS=$2
fi

# 15/02/22 - SJDK - Added the swif2 workflow as a variable you can specify here
#Workflow="cafe_aug08_${USER}" # Change this as desired
Workflow="cafe_replay_${USER}" # Change this as desired

# Input run numbers, this just points to a file which is a list of run numbers, one number per line
inputFile="${inputFile}${RunList}"

while true; do
    read -p "Do you wish to begin a new batch submission? (Please answer yes or no) " yn
    case $yn in
        [Yy]* )
            i=-1
            (
                # Finds number of lines of input file##   
                numlines=$(eval "wc -l < ${inputFile}")
		# Reads in input file, line by line
		while IFS='' read -r line || [[ -n "$line" ]]; do
                    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                    echo "Run number read from file: $line"
                    echo ""
                    # Each line is just a run number, set the variable to be the content of the line
                    runNum=$line
		    
		  
			
		    # Tape stub, you can point directly to a taped file and the farm job will do the jgetting for you, don't call it in your script! 
		    # The stub is chosen based upon run number. This is very rough, bit generally correct
		    if [ $runNum -ge 11599 ] && [ $runNum -le 16916 ]; then
			MSSstub='/mss/hallc/c-pionlt/raw/shms_all_%05d.dat'     # for cafe aug 08 data (taken under pion lt)
		    elif [ $runNum -ge 16917 ] && [ $runNum -le 17144 ]; then 
			MSSstub='/mss/hallc/c-cafe-2022/raw/shms_all_%05d.dat'  # cafe sep 18-29 2022 data
		    elif [ $runNum -ge 20776 ] && [ $runNum -le 20835 ]; then
			MSSstub='/mss/hallc/c-deuteron/raw/shms_all_%05d.dat'   # cafe feb 22-24 2023 data
		    elif [[ $runNum -lt 10000 ]]; then
			MSSstub='/mss/hallc/spring17/raw/coin_all_%05d.dat'
		    fi
		    # Output batch job text file, this is the script that is submitted as part of the job, 
		    # change the name of this as you want. Preferably, you should have different script name for each job so that you don't get any overwriting weirdness##
		    batch="${USER}_${runNum}_Job.txt"
		    tape_file=`printf $MSSstub $runNum`
		    # Print the size of the raw .dat file (converted to GB) to screen. sed command reads line 3 of the tape stub without the leading size=
		    TapeFileSize=$(($(sed -n '4 s/^[^=]*= *//p' < $tape_file)/1000000000))
		    if [[ $TapeFileSize == 0 ]];then
			TapeFileSize=2
		    fi
		    echo "Raw .dat file is "$TapeFileSize" GB"
		    echo "Job $(( $i + 2 ))/$(( $numlines +1 ))"
		    echo "Running ${batch} for ${runNum}"
		    cp /dev/null ${batch}
		    # Creation of batch script for submission, the script is just a series of commands
		    # We add these to the file (${batch}) via a series of piped echo commands
		    echo "PROJECT: c-comm2017" >> ${batch} # Or whatever your project is!
		    echo "TRACK: analysis" >> ${batch} ## Use this track for production running
		    #echo "TRACK: debug" >> ${batch} ### Use this track for testing, higher priority
		    echo "JOBNAME: cafe_${runNum}_${RunList}" >> ${batch} ## Change to be more specific if you want
		    # Request double the tape file size in space, for trunctuated replays edit down as needed
		    # Note, unless this is set typically replays will produce broken root files
		    echo "DISK_SPACE: "$(( $TapeFileSize * 2 ))" GB" >> ${batch}
		    if [[ $TapeFileSize -le 45 ]]; then # Assign memory based on size of tape file, should keep this as low as possible!
			echo "MEMORY: 3000 MB" >> ${batch}
		    elif [[ $TapeFileSize -ge 45 ]]; then
			echo "MEMORY: 4000 MB" >> ${batch}
		    elif [ $runNum -eq 16968  ]; then
			echo "MEMORY: 12000 MB" >> ${batch} 
		    fi
		    echo "CPU: 1" >> ${batch} ### hcana is single core, setting CPU higher will lower priority and gain you nothing!
		    echo "INPUT_FILES: ${tape_file}" >> ${batch}
		    echo "COMMAND:/w/hallc-scshelf2102/c-cafe-2022/cyero/cafe_offline_replay/replay_cafe_prod.sh ${runNum} ${MAXEVENTS}"  >> ${batch} ### Insert your script at end!
		    
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
		    i=$(( $i + 1 ))
		    if [ $i == $numlines ]; then
			echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			echo " "
			echo "###############################################################################################################"
			echo "############################################ END OF JOB SUBMISSIONS ###########################################"
			echo "###############################################################################################################"
			echo " "
		    fi
		    
		done < "$inputFile"
	    )
	    eval 'swif2 run ${Workflow}'
	    break;;
        [Nn]* ) 
	    exit;;
        * ) echo "Please answer yes or no.";;
    esac
done


