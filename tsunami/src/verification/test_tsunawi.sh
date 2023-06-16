#! /bin/bash

# Please check if these settings fit to your system,
# in particular the runcommand.

exe=ompTsuna.x
exepath=/home/users/ekrishnasamy/TsunAWI/tsunami/src
# runcommand="srun -c6 -t00:10:00 --qos=short -pmini,smp"   # Systems with Slurm, e.g. Ollie@AWI
runcommand="mylaptop"           # most systems, e.g., my Laptop 

# Defaults
# toggle this script, can be changed by comand line options
runref=0
comparison=0
verbose=0
runall=1

#########################################

usage()
{
    echo
    echo "Usage:  $(basename $0) [OPTIONS]"
    echo " "
    echo "  with the following OPTIONS:"
    echo "   -r      Run test cases and provide (new) reference results"
    echo "   -c      Run test cases and compare to reference results"
    echo "   -v      Be a bit verbose"
    echo "   -V<n>   Set verbosity to level <n>, 0=none, 1=a bit, 2=chatty"
    echo " "
    echo " All subdirectories with namelist.tsunami will be treated as test cases."
    echo " If you want to exclude tests, move the directory to another place."
    echo " Well, this is not an elegant solution and it should be possible"
    echo " to pass a list of tests to be performed..."
    echo " "

    exit 1
}
###

run_testcase()
{
    cd $1
    [[ $verbose -gt 0 ]] && echo "== Starting test case $exp"

    # Here and there, extra handling is needed. I have no idea how to make it more general...
    # Okushiri benchmark breaks if output at tide gauges is already present:
    if [[ "$exp" =~ "bm2_okushiri" ]]; then rm -f BM_oku_bm2_ch579.out; fi

    # Remove old checksums, if present. Otherwise, we might end up with failing
    # TsunAWI but do not recognize it because we work with old checksums.
    rm -f *_checksums.out

    # To test restarts, we have to run TsunAWI twice. And we must remove old restart
    # files from the last test.
    if ( grep -iq "write_restart=.true."  namelist.tsunami ); then
	rm -f RESTART*
	if [[ $verbose -gt 1 ]]; then
	    echo ""
	    echo "Starting experiment $exp"
	    $runcommand $exepath/$exe
	    echo "Restarting experiment $exp"
	    $runcommand $exepath/$exe
	else
	    $runcommand $exepath/$exe > /dev/null 2>&1
	    $runcommand $exepath/$exe > /dev/null 2>&1
	fi
    else
	if [[ $verbose -gt 1 ]]; then
	    echo ""
	    echo "Starting experiment $exp"
	    $runcommand $exepath/$exe
	else
	    $runcommand $exepath/$exe > /dev/null 2>&1
	fi
    fi
	

    # Ok, now we should have a new file *_checksums.out
    # The following line assumes that we have exactly one file named *_checksums.out!
    local checkfile=$(ls *_checksums.out)
    local reffile=${checkfile%.out}.ref
    
    # Is this a reference run?
    if [[ $runref -eq 1 ]]; then
	# If an old reference result is available, rename it
	[[ -f $reffile ]] && mv "$reffile" "${reffile%.ref}.bkp"
	# And now, rename output file to reference file
	mv "$checkfile" "$reffile"

    elif [[ $comparison -eq 1 ]]; then
	# This is a comparison run
	if [[ ! -f "$checkfile" ]]; then
	    echo " "
	    echo "!!! Test $exp failed - no checksums written !!!"
	    echo " "	
	elif [[ ! -f "$reffile" ]]; then
	    echo "!!!  No reference checksums for $exp !!!"
	else
	    compare "$checkfile" "$reffile"
	fi
    fi
    [[ $verbosity -gt 0 ]] && echo " "
    
    cd ..
}

###
compare()
{
    local checkfile=$1
    local reffile=$2

    if [[ $verbose -gt 0 ]]; then
	if ! diff -q "$checkfile" "$reffile" &> /dev/null; then
	    echo "$checkfile and $reffile differ"
	    echo "See $summary_file for details"
	else	    
	    echo "$checkfile and $reffile do not differ"
	fi
    fi

# As a start, just write new checksum, reference checksum and difference in one table
    echo "=== $exp ===" >> ../"$summary_file"
    echo -e "check \t\t  new value \t\t  reference value \t\t  difference \t\t relative (check-ref)/ref"  >>  ../"$summary_file"
    paste $checkfile $reffile | awk '{print $1, "\t",$2, "\t",$4, "\t",($2-$4),"\t",($2-$4)/$4}' >>  ../"$summary_file"

    echo " " >> ../"$summary_file"
    
    
}

##################################################################################

# And here we go. Check the arguments.

if [[ $# == 0 ]]; then usage; fi 


while getopts "chrvV:" OPT ; do
    case $OPT in
        c) comparison=1 ;;
	r) runref=1;  ;;
        v) verbose=1 ;;
        V) verbose=$OPTARG ;;
        h) usage ;;
        *) usage ;;
#
# one-liner by Malte to pass a list of scenarios (seperated by ":")
# 	 l) runall=( echo $OPTARG | awk -F: '{for(i=1;i<=NF;++i) print $i }'  ) ;;
# Later something like:
#  for ((i=0;i<${#runall[@]};++i)); do echo ${runall[$i]}; done

    esac
done
OPTIND=1  # reset for further usage

# Check, if any action is triggered
if [[ $comparison -eq 0 && $runref -eq 0 ]]; then
    echo " "
    echo "Please choose an action."
    usage
fi

# Check prerequisites
if [[ ! -d $exepath ]]; then
    echo "Sorry, $exepath is no directory"
    exit 1
fi
if [[ ! -x $exepath/$exe ]]; then
    echo "Sorry, $exepath/$exe is no executable"
    exit 1
fi
 
# Prepare for reference run
if [[ $runref -eq 1 ]]; then
# Check if older reference results shall be overwritten
    echo " "
    echo "Option -r is given. If reference checksums already exist, they will be renamed to <testid>_checksums.bkp."
    echo "Are you sure? (y/N)"; read YN
    if [[ $YN != "Y" &&  $YN != "y"  &&  $YN != "yes" ]]; then  exit 1; fi   
fi

# Prepare for comparison run
if [[ $comparison -eq 1 ]]; then
# Name the summary file with a unique name
    summary_file="summary_$(date '+%Y-%m-%d_%H-%M').out"
    if [[ -f $summary_file ]]; then          summary_file="summary_$(date '+%Y-%m-%d_%H-%M-%S').out"; fi
    if [[ -f $summary_file ]]; then sleep 1; summary_file="summary_$(date '+%Y-%m-%d_%H-%M-%S').out"; fi
    [[ $verbose -gt 0 ]] && echo "Output file: $summary_file"
fi

# Start the test cases   
if [[ $runall -eq 1 ]]; then
    # Cycle through all directories with test cases, i.e., all directories that contain a namelist.tsunami
    for exp in *; do
	if [[ -f $exp/namelist.tsunami ]]; then
	    run_testcase $exp
	fi
    done
fi

exit 0
