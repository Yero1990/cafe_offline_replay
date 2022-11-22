#! /bin/bash 


#filename='calib_test_runs.txt' 
#filename='h_singles_aug08.txt'
filename='optics_aug08.txt'
#filename='hms_xem.txt'

for run in $(cat $filename) ; do  
    
    #run=16036
    evt=200000
    
    # generic filename to be read by analzye_script
    root_file="ROOTfiles/cafe_replay_checkCalib_${run}_${evt}.root"
    
    # analys_script arguments 
    hms_pid="p"
    shms_pid="e"
    
    # define scripts
    replay_script="replay_cafe_checkCalib.sh"
    analyze_script="checkCalib.C"
    
    # define commands to run the scripts
    replay_cmd="./${replay_script} ${run} ${evt}"
    analyze_cmd="root -l -q -b \"${analyze_script}(\\\"${root_file}\\\", ${run}, \\\"${hms_pid}\\\", \\\"${shms_pid}\\\")\""  
    
    # ecevute the script
    echo "${replay_cmd}"
    eval ${replay_cmd}
    echo "${analyze_cmd}" 
    eval ${analyze_cmd}

done
