#!/bin/bash
#
# SConstruct phony targets run mu2e exe's via this script
# the one argument is the scons target name
#

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
    mkdir -p ${OUT}gen/txt
    scons $SC  -Q --tree=prune | deps -i > ${OUT}gen/txt/deps.txt
    [ $? -ne 0 ] && RC=1
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
