#! /bin/bash
#
# Check the status of packages of interest
#
# $Id: status.sh,v 1.3 2012/08/02 15:10:14 greenc Exp $
# $Date: 2012/08/02 15:10:14 $
# $Author: greenc $
#
# Contact person Rob Kutschke
#

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
for package in ${list}
do
  echo "Checking status of: " ${package}
  svn status BaBar/${package}
done
