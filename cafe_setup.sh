#!/bin/bash

# user input
fsys=$1

display_help()
{
    # Display Help
    echo "-------------------------------------------------------"
    echo "This shell script automatically sets up the necessary symbolic" 
    echo "links (or dir.) for the CaFe experiment raw data and output based on which "
    echo "machine (ifarm, cdaq, local) the user is at."
    echo ""
    echo "Syntax: ./cafe_setup.sh "
    echo ""
    echo "options:"
    echo ""
    echo "For users on IFARM ONLY: "
    echo "-f    ONLY use this option if you are running this shell script on ifarm. "
    echo "      This option selects filesystem in in which to read/write data " 
    echo "      from the CaFe experiment. "
    echo ""
    echo "      The optional arguments are: volatile, work, group"
    echo ""
    echo "      volatile (/volatile/hallc/c-cafe-2022/  4 TB high quota; 2 TB guarantee): "
    echo "              Files are NOT backed-up. Use for all large file output from analysis or simulation jobs. "
    echo "              When guarantee threshold is exceeded it is possible to have files auto-cleaned up (removed). "
    echo ""
    echo "      work (/work/hallc/c-cafe-2022/  1 TB  quota): "
    echo "             Files are NOT backed up. Good place for analysis software, database files, etc.  "
    echo "             These files should be in GitHub or have similar backups."
    echo ""
    echo "      group (/group/c-cafe-2022/  100 GB  quota): "
    echo "             Backed up regularly. Best place for analysis/replay scripts; software that is "
    echo "             being actively developed, etc. (But, still use GitHub!)"
    echo ""
    echo "      See https://hallcweb.jlab.org/wiki/index.php/CaFe_Disk_Space " 
    echo "      for detailed information on each of these filesystems."  
    echo "        "
    echo "examples: ./cafe_setup.sh volatile"
    echo "          ./cafe_setup.sh work" 
    echo "          ./cafe_setup.sh group"   
    echo "-------------------------------------------------------"    
}

if [ "$#" -eq 0 ]; then
    display_help
    exit 1
fi

#--- most recent online raw data is stored here (and then copied to tape shortly after)---
coda_raw="/net/cdaq/cdaql1data/coda/data/raw"
coda_raw_copiedtotape="/net/cdaq/cdaql1data/coda/data/raw.copiedtotape"



# where CaFe raw data output to be replayed will be stored (.dat files)
# but these are NOT directly accessible, one would have to look for them in cache.
# tape_raw_dir="/mss/hallc/c-cafe-2022/raw"
tape_raw_dir="/mss/hallc/c-pionlt/raw"

# tape volume for analysis output (simulation or replay output you want to keep long-term)
tape_analysis_out="/mss/hallc/c-cafe-2022/analysis" 


#--- define cache allocations (ONLY for raw and analyzed online data storage) ---
# cafe Sep data
cache_raw_dir_cafe="/cache/hallc/c-cafe-2022/raw/"
# cafe Aug data (optics, heep elastics, was written in pion lt raw dir)
cache_raw_dir_pionlt="/cache/hallc/c-pionlt/raw"
# cafe online analysis output would have been written here
cache_analysis_out="/cache/hallc/c-cafe-2022/analysis/"


set_hcana_link()
{
    if [[ -z $HCANALYZER ]]; then	
	echo ""
	echo "Environment variable: $HCANALYZER does NOT exist. "
	echo "Please make sure to do: source setup.sh(csh) in hcana first. " 
	echo "Trying to source environment:"
	echo "source ../hcana/setup.csh "
	source ../hcana/setup.csh
	cd
    else
	echo ""
	echo "Creating hcana symbolic link now  . . ."
	ln -sf $HCANALYZER"/hcana"
	#ln -sf $HCANALYZER"/libHallC.so"
	#ln -sf $HCANALYZER"/libHallC.so.0.90.0"
	echo ""
    fi    
}

# set link to cache cafe raw data
set_raw_link()
{
    mkdir CACHE_LINKS/
    cd CACHE_LINKS/

    unlink CACHE_LINKS/cache_cafe
    echo "Creating symlink to /cache/hallc/c-cafe-2022/raw/"
    ln -sf $cache_raw_dir_cafe cache_cafe
	
    unlink CACHE_LINKS/cache_pionlt
    echo "Creating symlink to /cache/hallc/c-pionlt/raw/"
    ln -sf $cache_raw_dir_pionlt cache_pionlt
    
}


# initialize machine flags to 0
# (depending on where this script gets called, it will turn ON one of these)
ifarm_flg=0
cdaq_flg=0



if echo $HOSTNAME | grep -q "ifarm"; then
    ifarm_flg=1
elif echo $HOSTNAME | grep -q "cdaq"; then
    cdaq_flg=1
fi


#if [[ ifarm_flg==0 &&  cdaq_flg==0 ]]; then
#    echo "***************************************"
#    echo " Did not recognize remote machine. "
#    echo " Please run: ./cafe_setup.sh -help "
#    echo " for help in running this script."
#    echo "***************************************"
#fi



#=================================
# ifarm
# (off-line experiment analysis)
#
# =================================
if [[ ifarm_flg -eq 1 ]]; then

  
    source setup.csh    
    echo ""
    echo "Checking if necessary directories or symlinks exist in remote machine: " ${USER}"@"${HOSTNAME}". . ."
    echo ""
    
    # setup the symbolic links to hcana
    set_hcana_link      

    # set raw data links
    set_raw_link
    
    if [[ $fsys == "volatile" ]]; then	     
	echo ""
	echo 'Setting up symbolic links to volatile filesystem on ifarm . . .'
	base_dir_voli="/volatile/hallc/c-cafe-2022/"	

	echo "Creating dir $base_dir_voli$USER . . ."
	mkdir $base_dir_voli$USER

	unlink REPORT_OUTPUT
	echo "Creating dir and symlink to $base_dir_voli$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_voli$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_voli$USER"/REPORT_OUTPUT"

	unlink CAFE_OUTPUT
	echo "Creating dir and symlink to $base_dir_voli$USER/CAFE_OUTPUT . . ."
	mkdir $base_dir_voli$USER"/CAFE_OUTPUT"
	ln -sf $base_dir_voli$USER"/CAFE_OUTPUT" 
	mkdir -p "CAFE_OUTPUT/ROOT"
	mkdir -p "CAFE_OUTPUT/REPORT"
	mkdir -p "CAFE_OUTPUT/PDF" 

	unlink ROOTfiles
	echo "Creating dir and symlink to $base_dir_voli$USER/ROOTfiles . . ."
	mkdir $base_dir_voli$USER"/ROOTfiles"
	ln -sf $base_dir_voli$USER"/ROOTfiles"
	echo ""
	
    elif [[ $fsys == "work" ]]; then	     
	echo ""
	echo 'Setting up symbolic links to work filesystem on ifarm . . .'
	base_dir_work="/work/hallc/c-cafe-2022/"

	echo "Creating dir $base_dir_work$USER . . ."
	mkdir $base_dir_work$USER

	unlink REPORT_OUTPUT
	echo "Creating dir and symlink to $base_dir_work$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_work$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_work$USER"/REPORT_OUTPUT"

	unlink CAFE_OUTPUT
	echo "Creating dir and symlink to $base_dir_work$USER/CAFE_OUTPUT . . ."	
	mkdir $base_dir_work$USER"/CAFE_OUTPUT" 
	ln -sf $base_dir_work$USER"/CAFE_OUTPUT"
	mkdir -p "CAFE_OUTPUT/ROOT"
	mkdir -p "CAFE_OUTPUT/REPORT"
	mkdir -p "CAFE_OUTPUT/PDF" 

	unlink ROOTfiles
	echo "Creating dir and symlink to $base_dir_work$USER/ROOTfiles . . ."
	mkdir $base_dir_work$USER"/ROOTfiles"
	ln -sf $base_dir_work$USER"/ROOTfiles"
	
	
    elif [[ $fsys == "group" ]]; then	     
	echo 'Setting up symbolic links to group filesystem on ifarm . . .'
	base_dir_group="/group/c-cafe-2022/"

	echo "Creating dir $base_dir_group$USER . . ."
	mkdir $base_dir_group$USER

	unlink REPORT_OUTPUT
	echo "Creating dir and symlink to $base_dir_group$USER/REPORT_OUTPUT . . ."
	mkdir $base_dir_group$USER"/REPORT_OUTPUT"
	ln -sf $base_dir_group$USER"/REPORT_OUTPUT"

	unlink CAFE_OUTPUT
	echo "Creating dir and symlink to $base_dir_group$USER/CAFE_OUTPUT . . ."
	mkdir $base_dir_group$USER"/CAFE_OUTPUT"
	ln -sf $base_dir_group$USER"/CAFE_OUTPUT"
	mkdir -p "CAFE_OUTPUT/ROOT"
	mkdir -p "CAFE_OUTPUT/REPORT"
	mkdir -p "CAFE_OUTPUT/PDF" 

	unlink ROOTfiles
	echo "Creating dir and symlink to $base_dir_group$USER/ROOTfiles . . ."
	mkdir $base_dir_group$USER"/ROOTfiles"
	ln -sf $base_dir_group$USER"/ROOTfiles"


	echo ""

    fi
fi


#===============================
# cdaq cluster
# (online experiment analysis)
#===============================

if [[ cdaq_flg -eq 1 ]]; then

    echo "Checking if necessary directories or symlinks exist in remote machine: " ${USER}"@"${HOSTNAME}". . ."

    # source cafe_online_replay
    source setup.csh
    
    # setup the symbolic links to hcana
    set_hcana_link
    
    base_dir_cdaq="/net/cdaq/cdaql1data/cdaq/hallc-online-cafe2022"

    echo "Creating symlink to ${coda_raw}"
    ln -sf $coda_raw

    echo "Creating symlink to ${coda_raw_copiedtotape}"
    ln -sf $coda_raw_copiedtotape 

    echo "Creating symlink to $cache_raw_dir_cafe"
    ln -sf $cache_raw_dir_cafe cache_cafe
    
    echo "Creating symlink to $cache_raw_dir_pionlt"
    ln -sf $cache_raw_dir_cafe cache_pionlt

    echo "Creating dir and symlink to $base_dir_cdaq/REPORT_OUTPUT . . ."
    mkdir $base_dir_cdaq"/REPORT_OUTPUT"
    ln -sf $base_dir_cdaq"/REPORT_OUTPUT"
    
    echo "Creating dir and symlink to $base_dir_cdaq/ROOTfiles . . ."
    mkdir $base_dir_cdaq"/ROOTfiles"
    ln -sf $base_dir_cdaq"/ROOTfiles"

    echo "Creating dir and symlink to $base_dir_cdaq/HISTOGRAMS . . ."
    mkdir $base_dir_cdaq"/HISTOGRAMS"
    ln -sf $base_dir_cdaq"/HISTOGRAMS"

    echo "Creating dir and symlink to $base_dir_cdaq/CAFE_OUTPUT . . ."    
    mkdir $base_dir_cdaq"/CAFE_OUTPUT"
    ln -sf $base_dir_cdaq"/CAFE_OUTPUT"
    mkdir -p "CAFE_OUTPUT/ROOT"
    mkdir -p "CAFE_OUTPUT/REPORT"
    mkdir -p "CAFE_OUTPUT/PDF"
	
fi


#=============================
# local
# (the user local computer)
#=============================

# assume user is local if NOT on cdaq or ifarm
if [[ ifarm_flg==0 && cdaq_flg==0 ]]; then
    
    # source cafe_online_replay (usually shell script on local machine)
    source setup.sh
    
    # This function checks if necessary dir. exists, else it creates them 
    dir_arr=("raw" "ROOTfiles" "REPORT_OUTPUT" "HISTOGRAMS" "CAFE_OUTPUT", "CACHE_LINKS")
    
    echo "Checking if necessary directories or symlinks exist in local machine: " ${USER}"@"${HOSTNAME}". . ."
    
    # setup the symbolic links to hcana
    set_hcana_link
    
    
    for i in "${dir_arr[@]}"
    do     
	if [[ -L "$i" && -d "$i" ]]; then
	    cmd="ls -l $i"
	    echo "$i is a symlink to a directory and it exists:"
	    eval $cmd 
	elif [[ -d "$i" ]]; then
	    echo "/$i directory exists"	
	else
	    echo "$i symlink is broken or /$i dir does not exist. Creating $i directory now . . ."
	    
	    cmd="mkdir $i"
	    echo $cmd
	    eval $cmd
	    echo "done!"
	fi    
    done
fi


