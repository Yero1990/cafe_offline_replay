#! /bin/bash

filename='ca48_runlist.txt'


for run in $(cat $filename) ; do    
    echo ${run}
    script="./run_cafe_sample.sh ${run} bcm_calib"
    echo "Executing command: ${script}"     
    eval ${script}    

done

  
