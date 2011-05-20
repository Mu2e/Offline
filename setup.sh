#
# $Id: setup.sh,v 1.25 2011/05/20 22:19:29 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/20 22:19:29 $
#
# Original author Rob Kutschke
#
# Setup the environment to build a full release of the Mu2e software.
# This presumes that you have already established the Mu2e environment
# for the machine on which you are running.
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

# Define the directory in which this file lives as the root of a release.
export MU2E_BASE_RELEASE=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && echo \$PWD`
echo "Base release directory is: " $MU2E_BASE_RELEASE

# These are needed by FileInPath.
unset  FW_BASE
export FW_RELEASE_BASE=$MU2E_BASE_RELEASE
export FW_SEARCH_PATH=$FW_RELEASE_BASE/:$FW_DATA_PATH/

# Setup the framework and its dependent products
setup art v0_07_02 -qa2:debug

# Geant4 and its cross-section files.
setup geant4 v4_9_4_p01 -qgcc45

setup g4neutron v3_14
setup g4emlow v6_19
setup g4photon v2_1
setup g4radiative v3_3
setup g4abla v3_0

# Other libraries we need. 
setup heppdt v3_04_01 -qgcc45

# The build system.
setup scons v1_3_0a -qgcc45

# Search path for fcl files
export FHICL_FILE_PATH=.:fcl;

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

# For now this must come after setup_mu2e_project - need to refactor the FW_ stuff.  It is no longer used.
export MU2E_SEARCH_PATH=$FW_SEARCH_PATH

# Build the symlink directories for the BaBar code.
# Only do so if the BaBar package is checked out locally.
if [  -f "./BaBar/makeInclude.sh" ]
then
  ./BaBar/makeInclude.sh
fi
