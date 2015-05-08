#! /bin/bash

usage() {
echo "
   samNoChildren [OPTIONS]  PARENT_DATASET CHILD_DATASET
      -c print only file counts of input and output datasets
      -s print only sam file names, not the full /pnfs names
      -f FILE File contains a list of the sam names of the parent 
                   dataset to operate on, instead of the whole dataset
      -h print help

   A file may have pointers to parent files in a parent dataset.
   This script reports the members of the parent dataset that have
   no members of child dataset pointing to them. This can be useful for 
   finding what child processing is not done yet. Finding only the sam
   names should take a few seconds, adding the full /pnfs path can add
   minutes per thousand output files.

   Examples:
   samNoChildren -c cnf.mu2e.cd3-beam-g4s2.0505a.fcl log.mu2e.cd3-beam-g4s2.0505a.log
   samNoChildren  cnf.mu2e.cd3-beam-g4s2.0505a.fcl log.mu2e.cd3-beam-g4s2.0505a.log
   samNoChildren -f temp.txt log.mu2e.cd3-beam-g4s2.0505a.log
"
}

export COUNTS=""
export SAMONLY=""
export INPUT=""

while getopts csf:h OPT; do
    case $OPT in
        c)
            export COUNTS=yes
            ;;
        s)
            export SAMONLY=yes
            ;;
        f)
            export INPUT=$OPTARG
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            echo unknown option, exiting
            exit 1
            ;;
     esac
done

shift $((OPTIND-1))
if [ "$INPUT" != ""  ]; then
  PDS=FILE
  CDS=$1
else
  PDS=$1
  CDS=$2
fi

if [[ "$PDS" == "" || "$CDS" == "" ]]; then
  echo "ERROR - parent or child datasets not provided"
  usage
  exit 1
fi

export SAM_EXPERIMENT=mu2e

# will need sam_web_client, setup if not already there
[ -z "$SETUP_SAM_WEB_CLIENT" ] && setup sam_web_client 
if [ "`which samweb`" == "" ]; then
  echo "ERROR - samweb not found and could not be setup"
  exit 3
fi

export FPU=`mktemp`
if [ "$INPUT" != "" ]; then
  while read FF
  do
    TT=`samweb file-lineage children $FF | \
     awk -F. -vd=$CDS '{if($1"."$2"."$3"."$4"."$6==d) print $0}'`
    if [ "$TT" = "" ]; then
      echo $FF >> $FPU
    fi
  done < $INPUT

else
  samweb list-files \
   "dh.dataset=$PDS and not  isparentof:(dh.dataset=$CDS)" > $FPU
fi

if [ "$COUNTS" ]; then
  if [ "$INPUT" != "" ]; then
    NPDS=`cat $INPUT | wc -l`
  else
    NPDS=`samweb list-files --summary "dh.dataset=$PDS" | grep "File count" | awk '{print $3}'`
  fi
  NCDS=`samweb list-files --summary "dh.dataset=$CDS" | grep "File count" | awk '{print $3}'`
  NUPS=`cat $FPU | wc -l`
  echo "File Counts: parent: $NPDS   child: $NCDS   parentNoChild: $NUPS"
  rm -f $FPU
  exit 0
fi

if [ "$SAMONLY" ]; then
  cat $FPU
  rm -f $FPU
  exit 0
fi

export FPUF=`mktemp`
N=0
while read FN
do
  FS=`samweb locate-file $FN | grep enstore | tr ":(" "  " | awk '{print $2}'`
  echo $FS/$FN >> $FPUF

#  N=$(($N+1))
#  T=$(($N % 100))
#  if [ $T -eq 0 ]; then
#    echo "Processing $N"
#  fi

done < $FPU

cat $FPUF
rm -f $FPU $FPUF
exit 0

#cnf.mu2e.cd3-beam-g4s2.0505a.fcl dataset that do not have children in the
#log.mu2e.cd3-beam-g4s2.0505a.log dataset.
