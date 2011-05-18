#! /bin/sh
#
# Populate the sym link directory neeed to use the BaBar style includes.
# Execute from Offline top level directory.
#

cd BaBar
if [ -e include ]
 then
  cd ..
  exit
fi

mkdir include
cd include

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for file in ${list}
do
  echo "checking ${file}"
  if [ ! -e ${file} ]
     then
     ln -s ../${file}/include ${file}
     echo "creating ${file}"
  fi
done

Dchlist="DchCalib DchData DchGeom DchGeomBase"
for file in ${Dchlist}
do
  echo "checking ${file}"
  if [ ! -e ${file} ]
    then
     ln -s ../Dch/${file}/include ${file}
     echo "creating ${file}"
  fi
done

# Two special cases.
if [ ! -e ErrLogger ]
  then
  ln -s  ../BaBar/include ErrLogger
fi
if [ ! -e PDT ]
  then
  ln -s  ../BaBar/include PDT
fi

cd ..
