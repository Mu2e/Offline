#
# $Id:
# $Author: wb $
# $Date: 2011/05/18 02:27:20 $
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
if [ "${FW_RELEASE_BASE}" = '' ];then
    echo "The environment variable FW_RELEASE_BASE is not set."
    echo "  It appears that you have not setup a base release of the Mu2e Offline software."
    echo "  You need to:"
    echo "     source /grid/fermiapp/mu2e/Offline/vx_y_z/setup.sh"
    echo "  where vx_y_z is the name of a tagged release."
    exit
fi

# Check that we are not stepping on an existing local setup.
if [ "${FW_BASE}" != '' ];then
    echo "The environment variable FW_BASE is already set."
    echo "So you already have a local environment established."
    echo "If you really want to do this, unset FWBASE and rerun this script."
    exit
fi

# Do the real work.
# Add the local directory to the search path.
# The / before the : in FW_SEARCH_PATH is significant.
export FW_BASE=$PWD
export FW_SEARCH_PATH=$PWD/:$FW_SEARCH_PATH
export MU2E_TEST_RELEASE=$PWD
