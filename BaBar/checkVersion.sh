#! /bin/bash
#
# Check the version of the checked out BaBar svn code and compare
# it to the required version, which is passed in as an argument.
#
# The check is done for each top level directory, not for each file.
# I believe that this is good enough.
#
# The check is done in two steps:
#  1) Check that the version of the checked out code matches the required version.
#  2) Check that there have been no changes since checkout.
#
# This code only prints warnings; it does not do an update.
#
# $Id: checkVersion.sh,v 1.5 2014/01/06 21:06:55 kutschke Exp $
# $Date: 2014/01/06 21:06:55 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

# The required version number is given as an input argument.
requiredVersion=$1

# If any second argument is present it enables verbose mode.
verbose=""
if [ ${#} = 2 ]; then
   verbose="verbose"
fi

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

list="BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra"
badBaBarPackage=""
modifiedPackages=""

if [ "${verbose}" != ''  ]; then
  echo "Checking for correct versions of BaBar packages."
fi

for package in ${list}
do

  # Step 1: checked out version matches the required version
  version=`svn info BaBar/${package} | awk '$1 == "Revision:" {print $2}'`

  if [ ${version} != ${requiredVersion} ]; then
    echo " "
    echo "ERROR " ${package} " does not match required version:"
    echo "   Required version: " ${requiredVersion}
    echo "   Actual version:   " ${version}
    echo " "
    set badBaBarPackage=${package} "  "  ${badBaBarPackage}
  fi

  # Step 2: no changes since checkout.
  count=`svn diff BaBar/${package} | wc | awk '{print $1}'`

  if [ "${count}" != "0" ]; then
    modifiedPackages=${modifiedPackages}"_"${package}
  fi

done

# Summary of step 1
if [ -z "${badBaBarPackage}" ]; then
  if [ "${verbose}" != ''  ]; then
    echo "All BaBar packages have the correct checked out version."
  fi
else
  echo " "
  echo "ERROR: Some BaBar packages have an incorrect version."
  echo "You must fix this before continuuing."
  echo "To update to the correct version issue the following command: "
  echo " "
  echo "BaBar/update.sh " ${requiredVersion}
  echo " "
fi

# Summary of step 2
if [ -z "${modifiedPackages}" ]; then
  if [ "${verbose}" != ''  ]; then
    echo "There have been no changes to the BaBar code since checkout."
  fi
else
  echo " "
  echo "WARNING:  Some BaBar packages have had changes since checkout."
  echo "          The packages are: " `echo ${modifiedPackages} | sed 's/\_/ /g'`
  echo "          If this is expected, carry on.  If not, contact an expert."
  echo " "
fi
