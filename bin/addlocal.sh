#
# $Id:
# $Author: gandr $
# $Date: 2012/09/24 19:00:42 $
#
# Original author Rob Kutschke
#
# Setup the environment to tell mute to look for certain files first in
# the local directory and then in the base release. This script requires
# that you have already setup a base release.  The files that can be looked
# in this manner include the geometry file(s), the generator configuration
# file(s), the conditions data file(s), the particle data tables and
# the magnetic field maps.  This list may be extended in the future.
#
# The code that opens input event-data files does not use this mechanism.
#

# Check that this is being used correctly
if [ "`basename $0 2>/dev/null`" = "localsetup.sh" ];then
    echo "You should be sourcing this file, not executing it."
    exit
fi

# Check that a base release has already been setup
if [ "${MU2E_BASE_RELEASE}" = '' ];then
    echo "The environment variable MU2E_BASE_RELEASE is not set."
    echo "  It appears that you have not setup a base release of the Mu2e Offline software."
    echo "  You need to:"
    echo "     source /grid/fermiapp/mu2e/Offline/vx_y_z/setup.sh"
    echo "  where vx_y_z is the name of a tagged release."
    return 1
fi

# Check that we are not stepping on an existing local setup.
if [ "${MU2E_BASE}" != '' ];then
    echo "The environment variable MU2E_BASE is already set."
    echo "So you already have a local environment established."
    echo "If you really want to do this, unset MU2E_BASEE and rerun this script."
    return 1
fi

# Do the real work.
# Add the local directory to the search path.
# The / before the : in FW_SEARCH_PATH is significant.
export MU2E_SATELLITE_RELEASE=$(/bin/pwd)
export MU2E_SEARCH_PATH=$(/bin/pwd)/:$MU2E_SEARCH_PATH
echo "MU2E_SEARCH_PATH:  " $MU2E_SEARCH_PATH

export FHICL_FILE_PATH=${MU2E_SATELLITE_RELEASE}:${MU2E_SATELLITE_RELEASE}/fcl:${FHICL_FILE_PATH}
