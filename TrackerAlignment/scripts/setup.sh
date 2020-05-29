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

setup millepede

python -m pip install --user -r ${TRKALIGN_SCRIPTS_DIR}/requirements.txt

# set up some convenience commands 

alias aligntrack_display='python ${TRKALIGN_SCRIPTS_DIR}/aligntrack_display.py ' 

function mu2ealign_genjobfcl() {
    cat job.fcl << EOF
#include "TrackerAlignment/fcl/cosmicAlign_timefit.fcl"
#include "JobConfig/reco/misalign_epilog.fcl"

services.ProditionsService.alignedTracker.useDb: true
services.ProditionsService.alignedTracker.verbose: 2

services.ProditionsService.mu2eDetector.useDb: true
services.ProditionsService.mu2eDetector.verbose: 2
physics.analyzers.TrkAnaNeg.diagLevel: 2

services.DbService.textFile: ["alignconstants_in.txt"]

physics.analyzers.AlignTrackCollector.diagLevel : 5

physics.analyzers.AlignTrackCollector.PlaneFilter : false
physics.analyzers.AlignTrackCollector.PlaneFilterList : [  ]

physics.analyzers.AlignTrackCollector.MilleFile : "MilleData.bin.gz"
physics.analyzers.AlignTrackCollector.GzipCompression : true 

physics.analyzers.AlignTrackCollector.SteerFile : "mp-steer.txt"
physics.analyzers.AlignTrackCollector.ParamFile : "mp-params.txt"
physics.analyzers.AlignTrackCollector.ConstrFile : "mp-constr.txt"

# diagnostics and plots, usually
services.TFileService.fileName: "TrackDiag.root"

EOF
    echo "Generated job.fcl!"

}

function mu2ealign() {
    COMMAND=$1

    if [[ $COMMAND == 'new']]; then
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
        mu2ealign_genjobfcl

        # produces a job.fcl to run and a seed alignment constant file
        # for DbService

        echo "Good to go!"
    fi
}

