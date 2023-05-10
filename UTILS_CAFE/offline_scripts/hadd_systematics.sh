#!/bin/bash

# script to combined multiple files for systematics study

#arguments pass to this bash script
#Use: on terminal type (for example) >> ./hadd_systematics Ca40 SRC

if [ $# -eq 0 ]; then
    echo ""
    echo "----------------------------------------------------------------------"
    echo "Use: ./hadd_systematics <target> <kin> "
    echo ""
    echo "<target>: Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197"
    echo "<kin>: MF, SRC"
    echo ""
    echo "description: reads a existing runlist: <target>_<kin>.txt 
          and generates an output runlist: <target>_<kin>_filelist.txt  
	  wich contain the actual list of  .root filenames 
	     filenames, which are used by hadd command to combine
             to a single .root file"
    echo "----------------------------------------------------------------------"
    echo ""

    exit 1
fi
    
#User input
tgt=$1   
kin=$2

HCREPLAY="/work/hallc/c-cafe-2022/$USER/cafe_offline_replay" 

# set generic path name of runlist to be read
filename="${HCREPLAY}/UTILS_CAFE/runlist/${tgt}_${kin}.txt"

# set output file containing the actual .root filenames to be combined
ofile="${tgt}_${kin}_filelist.txt"

# set generif output filename for the combined root files
combined_file="cafe_systematics_${tgt}_${kin}.root"

# make sure the output file does NOT exist, otherwise, delete and start over
if [ -e "$ofile" ]
then
    echo "file ${ofile} exists, will delete it and start over !"
    rm "$ofile"

    # loop over each run to create a list of runs
    for run in $(cat $filename) ; do    
	
	# set generic .root file name (of existing file)
	ifile="./cafe_prod_${tgt}_${kin}_${run}_-1_skimmed.root"

	# write the .root file names to file
	echo "$ifile" >> "$ofile"
      
    done

    # combine .root files
    haddCMD="hadd $combined_file @${ofile}"
    echo "Executing command: ${haddCMD}" 
    eval ${haddCMD}  
    
else
    echo "file ${ofile} does NOT exist, will create it !"

    # loop over each run to create a list of runs
    for run in $(cat $filename) ; do    
	
	# set generic .root file name (of existing file)  
	ifile="./cafe_prod_${tgt}_${kin}_${run}_-1_skimmed.root" 

	# write the .root file names to file
	echo "$ifile" >> "$ofile"
      
    done

    # combine .root files
    haddCMD="hadd $combined_file @${ofile}"
    echo "Executing command: ${haddCMD}" 
    eval ${haddCMD}  
    
fi



