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

# if running in muse, the default dir will be above Offline
if [ -n "$MUSE_WORK_DIR"  ]; then
    OUT=${MUSE_BUILD_BASE}/Offline/
    IN=Offline/
    SC=" -f $MUSE_DIR/python/SConstruct "
else
    OUT=""
    IN=""
    SC=""
fi

echo [$(date)] procs.sh starting $COMMAND

if [ "$COMMAND" == "DEPS"  ]; then
    # make a text file of the package dependencies
    if [ -n "$MUSE_WORK_DIR"  ]; then
	museDeps
	RC=$?
    else
	mkdir -p ${OUT}gen/txt
	scons $SC  -Q --tree=prune | deps -i > ${OUT}gen/txt/deps.txt
	[ $? -ne 0 ] && RC=1
    fi
elif [ "$COMMAND" == "GDML"  ]; then
    mkdir -p ${OUT}gen/gdml
    # if mu2e.gdml exists, save it
    if [ -f mu2e.gdml ]; then
	STR=$(date +%s)
	echo "save existing mu2e.gdml: mv mu2e.gdml mu2e.gdml.$STR"
	/bin/mv mu2e.gdml mu2e.gdml.$STR
    fi
    # make the standard gdml file
    mu2e -c ${IN}Mu2eG4/fcl/gdmldump.fcl; 
    [ $? -ne 0 ] && RC=1
    /bin/mv mu2e.gdml ${OUT}gen/gdml/mu2e.gdml
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "TEST03"  ]; then
    # see if this fcl runs
    mu2e -c ${IN}Mu2eG4/fcl/g4test_03.fcl -o /dev/null -T /dev/null
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "ROVERLAPS"  ]; then
    # fast root overlap check 
    TMP=$(mktemp)
    root -l -b -q ${IN}bin/overlapCheck.C\(\"${OUT}gen/gdml/mu2e.gdml\"\) >& $TMP
    [ $? -ne 0 ] && RC=1
    # a message on how many volumes checked
    grep "in Geometry imported from GDML" $TMP
    # check number of failures
    NI=$( grep "illegal" $TMP | awk '{print $NF}' )
    [ $NI -ne 0 ] && RC=1
    # also print them
    grep "illegal" $TMP
    cat $TMP | awk 'BEGIN{flag=0;}{if(flag==1) print $0; if($1=="===") flag=1; }'
    rm -f $TMP

elif [ "$COMMAND" == "GITPACK"  ]; then
    if [ -n "$MUSE_WORK_DIR"  ]; then
	cd Offline
    fi
    git repack -d -l
    [ $? -ne 0 ] && RC=1
    if [ -n "$MUSE_WORK_DIR"  ]; then
	cd ..
    fi
elif [ "$COMMAND" == "RMSO" ]; then
    if [ -n "$MUSE_WORK_DIR"  ]; then
	find ${MUSE_BUILD_BASE}/Offline -name "*.os" -delete
	find ${MUSE_BUILD_BASE}/Offline -name "*.o"  -delete
    else
	rm -rf tmp
	find . -name "*.os" -delete
	find . -name "*.o"  -delete
    fi
    rm -f .sconsign.dblite
elif [ "$COMMAND" == "VAL0" ]; then
    mkdir -p ${OUT}gen/val
    mu2e -n 5 -c ${IN}Validation/fcl/ceSimReco.fcl -o ceSimReco.art -T /dev/null
    [ $? -ne 0 ] && RC=1
    mu2e -s ceSimReco.art -T ${OUT}gen/val/ceSimReco_5000.root -c ${IN}Validation/fcl/val.fcl
    [ $? -ne 0 ] && RC=1
    rm -f ceSimReco.art
else
  echo "[$(date)] procs.sh did not parse argument: $@"
  exit 1
fi

if [ $RC -ne 0 ]; then
    echo "[$(date)] procs.sh failed $COMMAND"
    exit 1
fi

echo [$(date)] procs.sh success $COMMAND

exit 0
