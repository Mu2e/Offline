#
# $Id: setup.sh,v 1.57 2013/03/29 21:24:54 genser Exp $
# $Author: genser $
# $Date: 2013/03/29 21:24:54 $
#
# Original author Rob Kutschke
#
# Setup the environment to build or use a full release of the Mu2e software.
# This checks that you have already established the Mu2e environment.
#

if [ "`basename $0 2>/dev/null`" = "setup.sh" ];then
    echo "You should be sourcing this file, not executing it."
    exit 1
fi

if [ "${MU2E}" = '' ];then
    echo "The environment variable MU2E is not set."
    echo "You must setup the local Mu2e environment before sourcing this script."
    return 1
fi

# Protect against multiple invocation.
if [ "${MU2E_BASE_RELEASE}" != '' ];then
    echo "A base release has already been setup.  Hope that's OK."
    echo "The base release is: " ${MU2E_BASE_RELEASE}
    return 1
fi

# A very ill-defined state.  We have  a test release but no base release!
if [ "${MU2E_TEST_RELEASE}" != '' ];then
    echo "ERROR: A test release has already been setup but there is no base release."
    echo "Suggest that you log out, log in and restart from the beginning."
    echo "The test release is: " ${MU2E_TEST_RELEASE}
    return 1
fi

# Define the directory in which this file lives as the root of a release.
export MU2E_BASE_RELEASE=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && /bin/pwd`
echo "Base release directory is: " $MU2E_BASE_RELEASE

# Remove any test release environment.  TODO: test this and abort.
export MU2E_SEARCH_PATH=$MU2E_BASE_RELEASE/:$MU2E_DATA_PATH/
echo "MU2E_SEARCH_PATH:   "  $MU2E_SEARCH_PATH

# Setup the framework and its dependent products
setup art v1_03_08 -qe2:mu2e:prof

# Geant4 and its cross-section files.
setup geant4 v4_9_6_p01a -qe2:prof

# Other libraries we need.
setup heppdt v3_04_01 -qe2:prof
setup splines v1_01_03 -qe2:prof

# The build system.
setup scons v2_1_0 -qgcc47

# Search path for fcl files
export FHICL_FILE_PATH=${MU2E_BASE_RELEASE}:${MU2E_BASE_RELEASE}/fcl

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

# Check out the BaBar code.
# First build the symlink directory.  Then checkout the code.
if [  -f "./BaBar/makeInclude.sh" ]; then
  source ./BaBar/makeInclude.sh
  if [ ! -f "BaBar/BaBar/include/BaBar.hh" ]; then
   echo "Checking out the BaBar Kalman Filter code."
   ./BaBar/checkout.sh "-r 598"
  else
   echo "BaBar Kalman filter code already present. Not checking it out."
  fi
fi

# A hack that we hope can go away soon.
export G4LIBDIR=$G4LIB/$G4SYSTEM
