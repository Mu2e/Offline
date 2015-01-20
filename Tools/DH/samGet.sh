#! /bin/bash
#
# Ray Culbertson
#

usage() {
echo "
   samGet [OPTIONS] [-f FILE]  [-s FILEOFNAMES]  [-d DATASET]
      -n N    limit the nmber of files to retrieve
      -o DIR    direct output to directory DIR (deafult=.)
      -h print help

     Find certain files in SAM and copy them to a local directory.
     FILE is the comma-separated list of SAM names of file(s) 
      to be retrieved to the local directory
     FILEOFNAMES is a text file contains the sam names of files
     DATASET is the name of a dataset to retrieve. Since you probably don't 
     want all the files of a dataset, please limit the number with -n
     You need to "setup mu2e", kinit and getcert to run this procedure
     Only for interactive use, not intended for grid jobs.
"
}

export DIR=$PWD
export FILEOPT=""
export DSOPT=""
export FOFOPT=""
export NOPT=1000000
export OUTDIR=$PWD

while getopts d:f:n:o:s:h OPT; do
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
        n)
            export NOPT=$OPTARG
            ;;
        o)
            export OUTDIR=$OPTARG
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
# will need ifdh to do the copy
[ -z "$SETUP_IFDHC" ] && setup ifdhc

# write the sam file names to TMPF

TMPF=/tmp/${USER}_samgetf_$$
TMP=/tmp/${USER}_samget_$$
rm -f $TMPF $TMP

# for file names provided with -f
[ "$FILEOPT" != "" ] && echo $FILEOPT | tr "," "\n" >> $TMPF

# for file names provided with -s
[ "$FOFOPT" != "" ] && cat $FOFOPT >> $TMPF

# for file names provided with -d
if [ "$DSOPT" != ""  ]; 
then
  samweb list-files "dh.dataset=$DSOPT" >> $TMPF
fi

# translate sam file names to urls for ifdh

N=0
while read FN
do
echo samweb for file
  FNF=`samweb get-file-access-url --location=enstore $FN`
  echo $FNF " " $OUTDIR >> $TMP
  N=$(($N+1))
  if [ $N -ge $NOPT ]; then break; fi
done < $TMPF

ifdh cp -f $TMP
rm -f $TMP $TMPF

exit 0
