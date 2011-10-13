#! /bin/sh
#
# Checkout the packages of interest
#
# $Id: checkout.sh,v 1.1 2011/10/13 22:20:42 kutschke Exp $
# $Date: 2011/10/13 22:20:42 $
# $Author: kutschke $
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
