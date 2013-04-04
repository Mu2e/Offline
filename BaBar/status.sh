#! /bin/bash
#
# Check the status of packages of interest
#
# $Id: status.sh,v 1.4 2013/04/04 19:42:15 kutschke Exp $
# $Date: 2013/04/04 19:42:15 $
# $Author: kutschke $
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
  svn status -u BaBar/${package}
done
