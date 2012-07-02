#! /bin/bash
#
# Check the status of packages of interest
#
# $Id: status.sh,v 1.1 2012/07/02 15:44:08 kutschke Exp $
# $Date: 2012/07/02 15:44:08 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

cd BaBar

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for package in ${list}
do
  echo "Checking status of: " ${package}
  svn status ${package}
done

cd ..
