#! /bin/bash
#
# Ray Culbertson
#

usage() {
echo "
   samDatasets

     Print a list of mu2e datasets in SAM
"
}


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



# write the sam file names to TMPF

TMPF=/tmp/${USER}_samdsf_$$
TMP=/tmp/${USER}_samds_$$
rm -f $TMPF $TMP

samweb -e mu2e list-parameters dh.dataset | grep mu2e    >> $TMPF
samweb -e mu2e list-parameters dh.dataset | grep -v mu2e >> $TMPF

echo " files     MB      events     dataset"

while read DS
do

  rm -f $TMP
  samweb list-files --summary "dh.dataset=$DS" > $TMP
  NN=`cat $TMP | grep File | awk '{print $3}'`
  # in GB
  SZ=`cat $TMP | grep Total | awk '{printf "%d",$3/1000000.0}'`
  EC=`cat $TMP | grep Event | awk '{printf "%d",$3}'`

  #echo $DS NN $NN SZ $SZ  EC $EC
  if [ $NN -gt 0 ]; then
    printf "%6d %7d %10d   %s\n" $NN  $SZ $EC $DS
  fi

done < $TMPF

rm -f $TMP $TMPF


exit 0
