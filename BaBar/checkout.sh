#! /bin/sh
#
# Checkout the packages of interest
#
# $Id: checkout.sh,v 1.3 2012/09/26 18:14:45 brownd Exp $
# $Date: 2012/09/26 18:14:45 $
# $Author: brownd $
#
# Contact person Rob Kutschke
#

cd BaBar

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for package in ${list}
do
  svn co https://opteron05.lbl.gov/svn/MU2E/${package}/trunk ${package}
done

cd ..
