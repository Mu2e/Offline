#!/bin/bash
#
# Check the status of packages of interest
#
# $Id: update.sh,v 1.4 2012/08/02 15:10:14 greenc Exp $
# $Date: 2012/08/02 15:10:14 $
# $Author: greenc $
#
# Contact person Rob Kutschke
#

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

for PKG in BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra
do
  echo "Checking for updates in $PKG"
  svn update "BaBar/$PKG"
done
