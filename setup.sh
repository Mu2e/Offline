# Setup the environment to build or use a full release of the Mu2e software.
# This checks that you have already established the Mu2e environment.
#
# Original author Rob Kutschke
#

if [ "`basename $0 2>/dev/null`" = "setup.sh" ];then
    echo "You should be sourcing this file, not executing it." >&2
    exit 1
fi

if [[ ${#@} -ne 0 ]]; then
    cat >&2 <<EOF
ERROR: the setup.sh script does not accept arguments. Use

    $(dirname ${BASH_SOURCE})/buildopts

to query and define build configuraton.
EOF
    return 1
fi

if [ "${MU2E}" = '' ];then
    echo "The environment variable MU2E is not set." >&2
    echo "You must setup the local Mu2e environment before sourcing this script." >&2
    return 1
fi

# Protect against multiple invocation.
if [ "${MU2E_BASE_RELEASE}" != '' ];then
    echo "A base release has already been setup.  Hope that's OK." >&2
    echo "The base release is: " ${MU2E_BASE_RELEASE} >&2
    return 1
fi

# A very ill-defined state.  We have  a satellite release but no base release!
if [ "${MU2E_SATELLITE_RELEASE}" != '' ];then
    echo "ERROR: A satellite release has already been setup but there is no base release." >&2
    echo "Suggest that you log out, log in and restart from the beginning." >&2
    echo "The satellite release is: " ${MU2E_SATELLITE_RELEASE} >&2
    return 1
fi

# Define the directory in which this file lives as the root of a release.
export MU2E_BASE_RELEASE=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && /bin/pwd`
echo "Base release directory is: " $MU2E_BASE_RELEASE

# Remove any satellite release environment.  TODO: test this and abort.
export MU2E_SEARCH_PATH=$MU2E_BASE_RELEASE/:$MU2E_DATA_PATH/
echo "MU2E_SEARCH_PATH:   "  $MU2E_SEARCH_PATH

# Export a description of the configuration that we set up so that the
# build system can check consistency between setup and build configs.
export MU2E_SETUP_BUILDOPTS="$($MU2E_BASE_RELEASE/buildopts)"

build=$($MU2E_BASE_RELEASE/buildopts --build)

# This is the string to be used with the ups setup command for
# products that need qualifiers.  Note it includes the '+' character
# and is therefore different from the value shown in
# SETUP_<productname> environment vars, or by the "ups active" command.
export MU2E_UPS_QUALIFIERS=+e10:+${build}
export MU2E_ART_SQUALIFIER=s41

MU2E_G4_GRAPHICS_QUALIFIER=''
if [[ $($MU2E_BASE_RELEASE/buildopts --g4vis) == qt ]]; then
    MU2E_G4_GRAPHICS_QUALIFIER=':+qt'
fi

# Setup the framework and its dependent products
setup -B art v2_06_02a -q${MU2E_UPS_QUALIFIERS}

# root6 needs a path to include files to prevent some runtime warnings
export ROOT_INCLUDE_PATH=`dropit -s -p$ROOT_INCLUDE_PATH $MU2E_BASE_RELEASE`

# The interface to SAM - conflicts with ifdhc from the grid runtime environment
#setup -B ifdh_art v1_6_0 -q+e6:+${build}:+s5

# Geant4 and its cross-section files.
if [[ $($MU2E_BASE_RELEASE/buildopts --trigger) == "off" ]]; then
  setup -B geant4 v4_10_2_p03 -q${MU2E_UPS_QUALIFIERS}${MU2E_G4_GRAPHICS_QUALIFIER}
fi

# Get access to raw data formats.
# For now this is not generally available; once it is we can remove the conditional.
if [[ $($MU2E_BASE_RELEASE/buildopts --trigger) == "on" ]]; then
  setup -B mu2e_artdaq_core v1_00_06 -q${MU2E_UPS_QUALIFIERS}:+${MU2E_ART_SQUALIFIER}
fi

# Other libraries we need.
setup -B heppdt v3_04_01f -q${MU2E_UPS_QUALIFIERS}
setup -B BTrk   v1_01_08  -q${MU2E_UPS_QUALIFIERS}

# The build system.
setup -B scons v2_5_1 -q p2713b

# The debugger
setup -B gdb v7_12

# Search path for fcl files
export FHICL_FILE_PATH=${MU2E_BASE_RELEASE}:${MU2E_BASE_RELEASE}/fcl

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

#
if [ "${MU2E_BASE_RELEASE}" != `/bin/pwd` ]; then
  source ${MU2E_BASE_RELEASE}/bin/addlocal.sh
fi

# Environment variables used by the test build system.
export PACKAGE_SOURCE=${MU2E_BASE_RELEASE}
export BUILD_BASE=${MU2E_BASE_RELEASE}
