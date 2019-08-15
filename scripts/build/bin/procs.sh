#!/bin/bash
#
# SConstruct phony targets run mu2e exe's via this script
# the one argument is the scons target name
#

COMMAND=$1
RC=0

echo [$(date)] procs.sh starting $COMMAND

if [ "$COMMAND" == "DEPS"  ]; then
    # make a text file of the package dependencies
    scons -Q --tree=prune | deps -i > gen/txt/deps.txt
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "GDML"  ]; then
    # make the standard gdml file
    mu2e -c Mu2eG4/fcl/gdmldump.fcl; 
    [ $? -ne 0 ] && RC=1
    mv mu2e.gdml gen/gdml/mu2e.gdml
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "TEST03"  ]; then
    # see if this fcl runs
    mu2e -c Mu2eG4/fcl/g4test_03.fcl -o /dev/null -T /dev/null
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "ROVERLAPS"  ]; then
    # fast root overlap check 
    TMP=$(mktemp)
    root -l -b -q bin/overlapCheck.C\(\"gen/gdml/mu2e.gdml\"\) >& $TMP
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
    # set the remote to the writeable url
    git remote set-url origin  ssh://p-mu2eofflinesoftwaremu2eoffline@cdcvs.fnal.gov/cvs/projects/mu2eofflinesoftwaremu2eoffline/Offline.git
    [ $? -ne 0 ] && RC=1
    # make sure .git is packed
    git repack -d -l
    [ $? -ne 0 ] && RC=1
elif [ "$COMMAND" == "RMSO" ]; then
    rm -rf tmp
    rm -f .sconsign.dblite
    find . -name "*.os" -delete
    find . -name "*.o"  -delete
elif [ "$COMMAND" == "VAL0" ]; then
    mu2e -n 5000 -c Validation/fcl/ceSimReco.fcl -o ceSimReco.art -T /dev/null
    [ $? -ne 0 ] && RC=1
    mu2e -s ceSimReco.art -T gen/val/ceSimReco_5000.root -c Validation/fcl/val.fcl
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
