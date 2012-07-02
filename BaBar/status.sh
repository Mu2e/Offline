#! /bin/bash
#
# Check the status of packages of interest
#
# $Id: status.sh,v 1.2 2012/07/02 15:50:50 kutschke Exp $
# $Date: 2012/07/02 15:50:50 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for package in ${list}
do
  echo "Checking status of: " ${package}
  svn status BaBar/${package}
done
