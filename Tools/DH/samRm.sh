#! /bin/bash
#
# Ray Culbertson
#

################################
# usage function
################################
usage() {
echo "
   samRm [OPTIONS] [-f FILE]  [-s FILEOFNAMES] [-d DATASET]
      -n interpret file lists, but don't actually do the delete
      -h print help
  
     Permanently remove files from SAM database and 
     retire them in enstore (the tape system).
     You can re-upload a new file of the same name.
     The space on tape may eventually be overwritten.
     FILE is a comma-separated list of SAM file names to be deleted
     FILEOFNAMES is a text file contains the sam names of files
     DATASET is the name of a dataset to delete
     You need to "setup mu2e", kinit, and getcert to run this procedure
"
}

################################
# file delete function
################################
rmFile() {

  FN=`basename $1`

  if ! LOCS=`samweb locate-file $FN` ; then
     echo "ERROR: could not run locate on $FN"
    exit 1  
  fi


  N=`echo $LOCS | wc -w`
  if [ $N -gt 1 ]; then
    echo "ERROR: unexpected multiple locations for $1
$LOCS"
    exit 1
  fi

  if [ $LOCS ]; then
    MEDIA=`echo $LOCS | awk -F: '{print $1}'`
    PNFS=` echo $LOCS | awk -F: '{print $2}' | sed -e 's/([^()]*)//g' `
    PNFS=$PNFS/$FN

    if [ "$MEDIA" != "enstore" ]; then
      echo "ERROR: unknown value of MEDIA for $1
MEDIA = $MEDIA
echo $LOCS"
      exit 1
    fi

    echo "deleting $PNFS"
    if ! rm -f $PNFS ; then
     echo "ERROR: error during rm $PNFS"
     exit 1
    fi
    
  fi # end has a location

  echo "retiring $FN"
  if ! samweb retire-file $FN ; then
   echo "ERROR: error during sam retire-file $FN"
   exit 1
  fi


  return 0
}


################################
### main
################################

export FILEOPT=""
export DSOPT=""
export FOFOPT=""
export DOIT="yes"

while getopts d:f:s:hn OPT; do
    case $OPT in
        f)
            export FILEOPT=$OPTARG
            ;;
        d)
            export DSOPT=$OPTARG
            ;;
        s)
            export FOFOPT=$OPTARG
            ;;
        h)
            usage
            exit 0
            ;;
        n)
            export DOIT="no"
            ;;
        *)
            echo unknown option, exiting
            exit 1
            ;;
     esac
done

# must be ready to use sam_web_client
# which needs a grid cert
if ! grid-proxy-info >& /dev/null ; then
  echo "ERROR - grid certificate not found.  Please:
  kinit
  getcert
  export X509_USER_CERT=/tmp/x509up_u\`id -u\`"
  exit 2
fi

export SAM_EXPERIMENT=mu2e

# will need sam_web_client, setup if not already there
[ -z "$SETUP_SAM_WEB_CLIENT" ] && setup sam_web_client 

TMP=/tmp/${USER}_samrm_$$
touch $TMP

# for file names provided with -f
[ "$FILEOPT" != "" ] && echo $FILEOPT | tr "," "\n" >> $TMP

# for file names provided with -s
[ "$FOFOPT" != "" ] && cat $FOFOPT >> $TMP

# for file names provided with -d
if [ "$DSOPT" != ""  ]; 
then
  samweb list-files "dh.dataset=$DSOPT" >> $TMP
fi

N=`cat $TMP | wc -l`

if [ $N -le 0 ]; then
  echo Exiting - No files found 
  usage
  rm -f $TMP
  exit 1
fi

echo "Will delete $N files in sam and enstore, examples:"
cat $TMP | head -5

if [ "$DOIT" != "yes" ]; 
then
  exit 0
fi


echo "Now pausing 10 seconds to give you the option to ctrl-c"
sleep 10

while read FN
do
  rmFile $FN
done < $TMP

rm -f $TMP

exit 0
