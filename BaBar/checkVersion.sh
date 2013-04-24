#! /bin/bash
#
# Check the version of the checked out BaBar svn code and compare
# it to the required version, which is passed in as an argument.
#
# The check is done on a single directory.  I believe that this is enough
# but I can't be certain that it won't change in the future.
#
# This code only prints warnings; it does not do an update.
#
# $Id: checkVersion.sh,v 1.3 2013/04/24 19:01:19 kutschke Exp $
# $Date: 2013/04/24 19:01:19 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

# The required version number is given as an input argument.
requiredVersion=$1

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
badBaBarPackage=""

echo "Checking for correct versions of BaBar packages. May take 30 seconds"

for package in ${list}
do

  # Holds the output of an svn status command.
  tmpFile="svntmp_"${package}
  if [ -e $tmpFile ]; then
    /bin/rm -f $tmpFile
  fi

  #echo "Checking status of: " ${package}
  svn diff -r ${requiredVersion} BaBar/${package} >& ${tmpFile}

  # Any output is an error
  nLines=`wc $tmpFile | awk '{print $1}'`

  if [ $nLines != 0 ]; then
    echo " "
    echo "ERROR " ${package} " does not match required version " ${requiredVersion}
    svn status -u -q BaBar/${package}
    echo " "
    badBaBarPackage=${package}" "${badBaBarPackage}
  else
    /bin/rm -f ${tmpFile}
  fi

done

if [ -z "${badBaBarPackage}" ]; then
  echo "All BaBar packages pass the version check."
fi
