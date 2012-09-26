#! /bin/sh
#
# Checkout the packages of interest
#
# $Id: checkout.sh,v 1.2 2012/09/26 18:12:50 brownd Exp $
# $Date: 2012/09/26 18:12:50 $
# $Author: brownd $
#
# Contact person Rob Kutschke
#

cd BaBar

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for package in ${list}
do
  svn co -r 564 https://opteron05.lbl.gov/svn/MU2E/${package}/trunk ${package}
done

cd ..
