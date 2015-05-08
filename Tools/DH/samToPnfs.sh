#! /bin/bash
#
# Ray Culbertson
#

usage() {
echo "
   samToPnfs DATASET

     List the full /pnfs filespec of all files in the dataset.
     Useful for grid jobs which need a SAM datatset, but use a file 
     list for input instead of SAM.
"
}

DS="$1"

if [ ! "$DS" ];
then
  echo ERROR- no input dataset
  echo
  usage
  exit 1
fi

NF=`echo $DS | awk -F. '{print NF}'`

if [ $NF -ne 5 ]; then
  echo ERROR- input dataset does not have 5 dot fields
  echo
  usage
  exit 1
fi

F1=`echo $DS | awk -F. '{print $1}'`
F2=`echo $DS | awk -F. '{print $2}'`
F3=`echo $DS | awk -F. '{print $3}'`
F4=`echo $DS | awk -F. '{print $4}'`
F5=`echo $DS | awk -F. '{print $5}'`
# sim.mu2e.example-beam-g4s1.1812a.art

find /pnfs/mu2e/*/$F1/$F2/$F3/$F4/*/* -name "${F1}.${F2}.${F3}.${F4}.*.${F5}"

exit 0
