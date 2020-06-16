#!/usr/bin/bash

if [[ $0 == $BASH_SOURCE ]]; then
    echo "Please source this script."
    exit 1
fi

if [ -z ${MU2E_BASE_RELEASE} ]; then
    echo "You must first set up a release."
    return 1
fi

export TRKALIGN_SCRIPTS_DIR="${MU2E_BASE_RELEASE}/TrackerAlignment/scripts"
export DS_COSMIC_NOFIELD_ALIGNSELECT="/pnfs/mu2e/persistent/users/mu2epro/MDC2020Dev/DS-cosmic-nofield-alignselect/sources.txt"

setup millepede

# python -m pip install --user -r ${TRKALIGN_SCRIPTS_DIR}/requirements.txt

# set up some convenience commands 

alias aligntrack_display='python ${TRKALIGN_SCRIPTS_DIR}/aligntrack_display.py ' 

function mu2ealign_genparallel() {
    rm job_part*.fcl > /dev/null

    nparts=$1
    nfiles=$2
    echo "generating fcls for $nparts processes and $nfiles total files using job.fcl"
    END=$nparts
    head -n $nfiles sources.txt > sources.txt.tmp
    for ((i=1;i<=END;i++)); do
        cp job.fcl job_part$i.fcl
        sed -i "s/MilleData.bin/MilleData.bin.${i}/g" job_part$i.fcl
        sed -i "s/mp-steer.txt/mp-steer.txt.${i}/g" job_part$i.fcl
        sed -i "s/mp-constr.txt/mp-constr.txt.${i}/g" job_part$i.fcl
        sed -i "s/mp-params.txt/mp-params.txt.${i}/g" job_part$i.fcl
        sed -i "s/TrackDiag.root/TrackDiag.root.${i}/g" job_part$i.fcl

        split --number=$i/$nparts -d sources.txt.tmp > sources_job_part${i}.txt
    done

    rm sources.txt.tmp
}

function mu2ealign_runjobs() {
    if [ ! -f "job_part1.fcl" ]; then
        mu2e -c job.fcl -S sources.txt > job.log 2>&1 &
        return 0
    fi

    i=0
    for f in job_part*.fcl; do
        i=$((i+1))
        mu2e -c job_part${i}.fcl -S sources_job_part$i.txt > job_part$i.log 2>&1  &
    done

    echo "Started $i jobs.. see job_part{job number}.log files for progress.."
}

function mu2ealign_checkcomplete() {
    if [ ! -f "job_part1.fcl" ]; then
        if ! grep -q "Art has completed and will exit with status 0." job.log ; then 
            echo "Job incomplete! See below (job.log):"
            echo "-------------------------------------------"
            tail -n 5 job.log
            echo "-------------------------------------------"
            return 1
        fi
        return 0
    fi

    i=0
    rc=0
    for f in job_part*.fcl; do
        i=$((i+1))
        if ! grep -q "Art has completed and will exit with status 0." job_part$i.log ; then 
            echo "Job $i incomplete! See below (job_part$i.log):"
            echo "-------------------------------------------"
            tail -n 5 job_part$i.log
            echo "-------------------------------------------"
            rc=1
        #else
        #    echo "Job $i complete!"
        #    grep "TimeReport" job_part$i.log
        #    echo ""
        fi
    done

    return $rc
}

function mu2ealign_mergeoutput() {
    if [ ! -f "job_part1.fcl" ]; then
        echo "nothing to merge"
        return 1
    fi

    if [ ! -f "mp-steer.txt" ]; then
        python ${MU2E_BASE_RELEASE}/TrackerAlignment/scripts/mergesteer.py mp-steer.txt.* > mp-steer.txt
    else
        echo "mp-steer.txt already exists - skip"
    fi

    if [ ! -f "TrackDiag.root" ]; then
        hadd -f TrackDiag.root TrackDiag.root.*
    else
        echo "TrackDiag.root already exists - skip"
    fi
}

function mu2ealign_genjobfcl() {
    cp ${MU2E_BASE_RELEASE}/TrackerAlignment/fcl/job_template.fcl job.fcl
    echo "Generated new job.fcl!"
    echo "Using DS_COSMIC_NOFIELD_ALIGNSELECT as dataset. ( 4 files )"
    echo "Please change sources.txt if you want to use something else."
    head -n 4 ${DS_COSMIC_NOFIELD_ALIGNSELECT} > sources.txt
}

function mu2ealign_progress() {
    tracks=0
    for f in job_part*.log; do
        tcount=$(tac $f | grep -m1 "wrote track" | awk 'NF>1{print $NF}' | tr -d '\n' | tr -d '\r')
        if [ "$tcount" = "" ]; then 
          tcount=0
        fi
        echo $tcount
        tracks=$((tracks + tcount))
    done
    echo -en "\r... $tracks tracks ..."
}

function mu2ealign_runNaligniters() {
    END=$1
    (
        echo "Working directory: $(pwd)"
        echo "Alignment track collection: iteration 0"

        # run first alignment iteration
        if [ ! -f "alignconstants_out.txt" ]; then
            mu2ealign run
            lastm2epid=$!
            while [ -d "/proc/$lastm2epid" ]; do
                mu2ealign_progress
                sleep 2
            done
            echo ""
            echo "... nearly done ..."

            wait;

            mu2ealign pede
            lastconsts=$(pwd)/alignconstants_out.txt
        else
            echo "Already complete - skip"
        fi

        for ((alignjobn=1;alignjobn<=END-1;alignjobn++)); do

            mkdir -p iter$alignjobn || return 1
            cd iter$alignjobn || return 1
            if [ ! -f "alignconstants_out.txt" ]; then
                echo "Working directory: $(pwd)"

                mu2ealign new $lastconsts

                echo "Alignment track collection: iteration $alignjobn"
                
                mu2ealign run 

                lastm2epid=$!
                while [ -d "/proc/$lastm2epid" ]; do
                    mu2ealign_progress
                    sleep 2
                done
                echo ""
                echo "... nearly done ..."
                wait 

                mu2ealign pede
            else
                echo "Already complete - skip"
            fi
            lastconsts=$(pwd)/alignconstants_out.txt
            cd ..
        done

        # run nominal config for comparison
        mkdir nominal
        cd nominal 
        echo "Working directory: $(pwd)"

        mu2ealign new ${MU2E_BASE_RELEASE}/TrackerAlignment/test/misalignments/nominal.txt

        rm *.fcl 
        rm sources*.txt 
        cp ../*.fcl .
        cp ../sources_*.txt .

        echo "Alignment track collection: NOMINAL geom"
            
        mu2ealign run 

        wait 
        mu2ealign_checkcomplete || return 1
        mu2ealign_mergeoutput


        echo "Complete! Final alignment constants are in $lastconsts"
    )
}

function mu2ealign() {
    COMMAND=$1

    if [[ $COMMAND == "new" ]]; then
        if [ "$(ls -A $PWD)" ]; then
            echo "Please re-run this command inside an empty directory."

            echo "i.e."
            echo "$ mkdir align_iter0 && cd align_iter0"
            echo "$ mu2ealign new <path to alignment constants file>"
            return 1
        fi

        if [ -z "$2" ]; then 
            echo "usage: "
            echo "$ mu2ealign new <alignment constants file>"
            return 1
        fi

        # generate a working directory in CWD
        # the job uses the alignment constants file in $2
        ALIGN_CONST_FILE=$2

        if [ ! -f "${ALIGN_CONST_FILE}" ]; then
            TESTFILE="${MU2E_BASE_RELEASE}/TrackerAlignment/test/misalignments/$2.txt"
            if [ -f "${TESTFILE}" ]; then 
                ALIGN_CONST_FILE=${TESTFILE}

                echo "using: ${ALIGN_CONST_FILE}"
            else
                echo "$ALIGN_CONST_FILE does not exist."
                return 1
            fi
        fi

        cp ${ALIGN_CONST_FILE} alignconstants_in.txt

        JOB_FCL_FILE=$(dirname ${ALIGN_CONST_FILE})/job.fcl

        if [ -f ${JOB_FCL_FILE} ]; then 
            # copy old fcl over
            cp $(dirname ${ALIGN_CONST_FILE})/*.fcl .
            cp $(dirname ${ALIGN_CONST_FILE})/sources*.txt .

            echo "Copied previous job configuration!"
        else
            mu2ealign_genjobfcl
        fi

        # produces a job.fcl to run and a seed alignment constant file
        # for DbService

        if [ ! -f "sources_job_part1.txt" ]; then
            echo "If you want to configure for multiple jobs, run mu2ealign parallel <NJOBS>"
        fi
        echo "Run 'mu2ealign run' to start. To run multiple alignment iterations, run mu2ealign autorun <NITERS>."

    elif [[ $COMMAND == "parallel" ]]; then 
        mu2ealign_genparallel $2 4

    elif [[ $COMMAND == "run" ]]; then 
        mu2ealign_runjobs

    elif [[ $COMMAND == "autorun" ]]; then
        mu2ealign_runNaligniters $2

    elif [[ $COMMAND == "pede" ]]; then
        # check completion of jobs
        mu2ealign_checkcomplete || return 1
        mu2ealign_mergeoutput

        pede mp-steer.txt || return 1

        if [ -f "millepede.res" ]; then 
            python ${TRKALIGN_SCRIPTS_DIR}/mp2prod.py > alignconstants_out.txt

            echo "Generated new alignment constants in alignconstants_out.txt."
        fi

    elif [[ $COMMAND == "help" ]]; then 
        echo "Available commands: "
        echo "mu2ealign new <path to alignment constants txt file>: create a new alignment working directory using specified alignment constants"
        echo "mu2ealign run: start the alignment track collection in the current working directory"
        echo "mu2ealign pede: check jobs have completed, merge results if needed and run PEDE. Produce new alignment constants"
    fi
}

