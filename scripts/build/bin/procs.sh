#!/bin/bash
#
# SConstruct phony targets run mu2e exe's via this script
# the one argument is the scons target name
#


#
# a function to make the deps file in the Muse case
#
museDeps() {

    local INCPATH=$( scons  -f $MUSE_DIR/python/SConstruct  -s INCPATH \
        | sed 's|'$PWD'|.|g' )

    mkdir -p $MUSE_BUILD_BASE/Offline/gen/txt
    local DFILE=$MUSE_BUILD_BASE/Offline/gen/txt/deps.txt
    rm -f $DFILE
    touch $DFILE

    # remove Project (unmaintained) and .git
    local DIRS=$( find Offline -mindepth 1 -maxdepth 1 -type d | sed 's|^./||' | \
        grep -v Projects | grep -v git )

    for DD in $DIRS
    do
        local REPO=$(echo $DD | awk -F/ '{print $2}')
        # remove CRVResponse wls subdir - wrong include paths
        local FILES=$(find $DD \( -name "classes.h" -or -name "*.hh" -or -name "*.cc"  \) | \
            grep -v wls )
        local N=$( echo $FILES | wc -w )
        if [ $N -gt 0 ]; then
            cat $FILES > sconsInclude.cc
            g++ --std=c++17  -M $INCPATH sconsInclude.cc \
                | sed -e 's|\\|\n|g' -e 's/ /\n/' | grep -v sconsInclude | \
                sed '/^\/usr\/*/d' | sed -e '/^[[:space:]]*$/d'   > sconsInclude.txt
            local RC=$?
            if [ $RC -ne 0 ]; then
                echo "ERROR compile error in $DD"
                rm -f sconsInclude*
                return 1
            fi
            local HDEPS=$(cat sconsInclude.txt | grep Offline | awk -F/ '{print $2}' | \
                sort | uniq | grep -v $REPO | tr "\n" " " )
            local PDEPS=$(cat sconsInclude.txt | grep cvmfs | awk -F/ '{print $5}' | \
                sort | uniq | grep -v $REPO | tr "\n" " ")
            local NI=$(cat sconsInclude.txt | wc -l)
            echo $REPO $N files $NI includes

            echo "HDR $REPO $HDEPS" >> $DFILE
            echo "PRD $REPO $PDEPS" >>  $DFILE
        else
            echo $REPO 0 files
        fi
    done

    rm -f sconsInclude*

    return 0
}



COMMAND=$1
RC=0

echo [$(date)] procs.sh starting $COMMAND

if [ "$COMMAND" == "DEPS"  ]; then
    # make a text file of the package dependencies
    museDeps
    RC=$?
elif [ "$COMMAND" == "GDML"  ]; then
    # find gdml fcl files and make all GDML files, add to build dir
    RC=0
    mkdir -p ${MUSE_BUILD_BASE}/Offline/gen/gdml
    FCLS=$(ls Offline/Mu2eG4/fcl/gdmldump_* )
    OWD=$PWD
    TWD=$(mktemp -d)
    cd $TWD
    for FCL in $FCLS
    do
        echo
        echo "Starting GDML for $FCL"
        mu2e -c $FCL
        # some tags and geometries are inconsistent
        # so do not fail on failure to make GDML
    done
    mv *.gdml ${MUSE_WORK_DIR}/${MUSE_BUILD_BASE}/Offline/gen/gdml
    cd $OWD
    rmdir $TWD
elif [ "$COMMAND" == "TEST03"  ]; then
    # see if this fcl runs
    mu2e -c Offline/Mu2eG4/fcl/g4test_03.fcl -o /dev/null -T /dev/null
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "ROVERLAPS"  ]; then
    # fast root overlap check
    Offline/bin/overlapCheck.sh -r -q -b
elif [ "$COMMAND" == "GITPACK"  ]; then
    git -C Offline repack -d -l
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "RMSO" ]; then
    find ${MUSE_BUILD_BASE}/Offline -name "*.os" -delete
    find ${MUSE_BUILD_BASE}/Offline -name "*.o"  -delete
    rm -f .sconsign.dblite
fi

if [ $RC -ne 0 ]; then
    echo "[$(date)] procs.sh failed $COMMAND"
    exit 1
fi

echo [$(date)] procs.sh success $COMMAND

exit 0
