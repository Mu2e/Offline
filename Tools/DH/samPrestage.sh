#! /bin/bash
#
# Ray Culbertson
#

usage() {
echo "
   samPrestage DATASETDEFINTION
     -v verbose
  
     Copy the files of a sam dataset defintion from tape to dCache.
     When the script completes, the files will be in the tape queue,
     not necessary in dCache yet.
     You need to "setup mu2e", kinit and getcert to run this procedure
     Only for interactive use, not intended for grid jobs.
"
}

export V=""
while getopts hv OPT; do
    case $OPT in
        h)
            usage
            exit 0
            ;;
        v)
            export V="true"
            shift;
            ;;
        *)
            echo unknown option, exiting
            exit 1
            ;;
     esac
done

if [ ! $1 ];
then
  usage
  exit 1
fi


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

export SAM_DD=$1

echo "Summary of files to be prestaged:"
samweb list-files --summary "dataset_def_name=$SAM_DD"

export SAM_PROJECT="prestage_${SAM_DD}_`date +%s`"
export SAM_PROJECT_URL=`samweb start-project --defname=$SAM_DD $SAM_PROJECT`

if [ $V ]; then
  echo SAM_PROJECT_URL=$SAM_PROJECT_URL
fi
if [ ! $SAM_PROJECT_URL ]; then
  echo SAM_PROJECT_URL not set - exiting!
  exit 1
fi

export SAM_CONSUMER_ID=`samweb start-process --appname=null --appversion=0 $SAM_PROJECT_URL`
if [ $V ]; then
  echo SAM_CONSUMER_ID=$SAM_CONSUMER_ID
fi
if [ ! $SAM_CONSUMER_ID ]; then
  echo SAM_CONSUMER_ID not set - exiting!
  exit 1
fi

N=0
export SAM_FILE=`samweb get-next-file $SAM_PROJECT_URL $SAM_CONSUMER_ID`
while [ $SAM_FILE ];
do 
  N=$(($N+1))
  if [ $V ]; then
    echo SAM_FILE=$SAM_FILE
  fi

  samweb release-file --status=ok  $SAM_PROJECT_URL $SAM_CONSUMER_ID $SAM_FILE

  TEST=$(($N%10))
  if [[ $N -lt 100 && $TEST -eq 0 ]]; then
    echo Processing $N
  fi
  TEST=$(($N%100))
  if [[ $N -lt 1000 && $TEST -eq 0 ]]; then
    echo Processing $N
  fi
  TEST=$(($N%1000))
  if [[ $TEST -eq 0 ]]; then
    echo Processing $N
  fi

  export SAM_FILE=`samweb get-next-file $SAM_PROJECT_URL $SAM_CONSUMER_ID`

done

echo "Prestaged $N files"

exit 0
