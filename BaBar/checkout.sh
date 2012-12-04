#! /bin/sh
#
# Checkout the packages of interest
#
# $Id: checkout.sh,v 1.5 2012/12/04 00:51:28 tassiell Exp $
# $Date: 2012/12/04 00:51:28 $
# $Author: tassiell $
#
# Contact person Rob Kutschke
#

revision=""
if [ $# == 1 ]; then
 revision=$1
 echo "Checking out svn revision: " ${revision}
fi

cd BaBar

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra Dch"
for package in ${list}
do
  svn co ${revision} https://opteron05.lbl.gov/svn/MU2E/${package}/trunk ${package}
done

cd ..
