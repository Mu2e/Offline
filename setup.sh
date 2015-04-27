#
# $Id: setup.sh,v 1.70 2014/08/22 17:45:41 brownd Exp $
# $Author: brownd $
# $Date: 2014/08/22 17:45:41 $
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

# A very ill-defined state.  We have  a satellite release but no base release!
if [ "${MU2E_SATELLITE_RELEASE}" != '' ];then
    echo "ERROR: A satellite release has already been setup but there is no base release."
    echo "Suggest that you log out, log in and restart from the beginning."
    echo "The satellite release is: " ${MU2E_SATELLITE_RELEASE}
    return 1
fi

# Define the directory in which this file lives as the root of a release.
export MU2E_BASE_RELEASE=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && /bin/pwd`
echo "Base release directory is: " $MU2E_BASE_RELEASE

# Remove any satellite release environment.  TODO: test this and abort.
export MU2E_SEARCH_PATH=$MU2E_BASE_RELEASE/:$MU2E_DATA_PATH/
echo "MU2E_SEARCH_PATH:   "  $MU2E_SEARCH_PATH

build=${1:-prof}
if [ "${build}" == "debug" ];then
    # echo "debug option selected; setting up gdb"
    setup gdb v7_8
    # echo "use scons --mu2elevel=debug to create debug libraries"
fi

# This is the string to be used with the ups setup command for
# products that need qualifiers.  Note it includes the '+' character
# and is therefore different from the value shown in
# SETUP_<productname> environment vars, or by the "ups active" command.
export MU2E_UPS_QUALIFIERS=+e7:+${build}

# Setup the framework and its dependent products
setup -B art v1_13_01 -q${MU2E_UPS_QUALIFIERS}

# The interface to SAM - conflicts with ifdhc from the grid runtime environment
#setup -B ifdh_art v1_6_0 -q+e6:+${build}:+s5

# Geant4 and its cross-section files.
setup -B geant4 v4_9_6_p04a -q${MU2E_UPS_QUALIFIERS}

# Other libraries we need.
setup -B heppdt v3_04_01c -q${MU2E_UPS_QUALIFIERS}
setup -B splines v1_06_00 -q${MU2E_UPS_QUALIFIERS}

# The build system.
setup -B scons v2_3_4


# Search path for fcl files
export FHICL_FILE_PATH=${MU2E_BASE_RELEASE}:${MU2E_BASE_RELEASE}/fcl

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

# Check out the BaBar code.
# First build the symlink directory.  Then checkout the code.
babarversion=616
if [  -f "./BaBar/makeInclude.sh" ]; then
  source ./BaBar/makeInclude.sh
  if [ ! -f "BaBar/BaBar/include/BaBar.hh" ]; then
   echo "Checking out the BaBar Kalman Filter code."
   ./BaBar/checkout.sh "-r ${babarversion}"
  else
    # Skip check during grid jobs or else they will DOS attack the svn repository.
    if [ -z "${PROCESS}" ] || [ "${PROCESS}" == 0 ]; then
      ./BaBar/checkVersion.sh ${babarversion}
    else
      echo "Grid job with PROCESS != 0 detected. Skipping version check of BaBar code."
    fi
  fi
fi

#
if [ "${MU2E_BASE_RELEASE}" != `/bin/pwd` ]; then
  source ${MU2E_BASE_RELEASE}/bin/addlocal.sh
fi

# A hack that we hope can go away soon.
export G4LIBDIR=$G4LIB/$G4SYSTEM

# Add useful functions to the shell environment.
source ${MU2E_BASE_RELEASE}/bin/functions.sh
