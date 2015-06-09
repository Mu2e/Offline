#!/bin/bash
#
# Check the status of packages of interest
#
# $Id: update.sh,v 1.6 2013/04/04 19:48:30 kutschke Exp $
# $Date: 2013/04/04 19:48:30 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

revision=""
if [ $# == 1 ]; then
 revision=$1
 echo "Updating to svn revision: " ${revision}
fi


for PKG in BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra Dch
do
  echo "Checking for updates in $PKG"
  svn update -r ${revision} "BaBar/$PKG"
done
