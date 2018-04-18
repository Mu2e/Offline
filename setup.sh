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

# Export a description of the configuration that we set up so that the
# build system can check consistency between setup and build configs.
export MU2E_SETUP_BUILDOPTS="$($MU2E_BASE_RELEASE/buildopts)"

build=$($MU2E_BASE_RELEASE/buildopts --build)

# This is the string to be used with the ups setup command for
# products that need qualifiers.  Note it includes the '+' character
# and is therefore different from the value shown in
# SETUP_<productname> environment vars, or by the "ups active" command.
export MU2E_UPS_QUALIFIERS=+e15:+${build}
export MU2E_ART_SQUALIFIER=s66

MU2E_G4_GRAPHICS_QUALIFIER=''
if [[ $($MU2E_BASE_RELEASE/buildopts --g4vis) == qt ]]; then
    MU2E_G4_GRAPHICS_QUALIFIER=':+qt'
fi

# Setup the framework and its dependent products
setup -B art v2_10_04 -q${MU2E_UPS_QUALIFIERS}

# Geant4 and its cross-section files.
if [[ $($MU2E_BASE_RELEASE/buildopts --trigger) == "off" ]]; then
  setup -B geant4 v4_10_2_p03e -q${MU2E_UPS_QUALIFIERS}${MU2E_G4_GRAPHICS_QUALIFIER}
else
  setup -B xerces_c v3_1_4b   -q${MU2E_UPS_QUALIFIERS}
fi

# Get access to raw data formats.
setup -B mu2e_artdaq_core v1_02_01f -q${MU2E_UPS_QUALIFIERS}:+${MU2E_ART_SQUALIFIER}:offline

# Other libraries we need.

setup -B heppdt   v3_04_01g -q${MU2E_UPS_QUALIFIERS}
setup -B BTrk   v1_02_11  -q${MU2E_UPS_QUALIFIERS}
setup -B cry   v1_7j  -q${MU2E_UPS_QUALIFIERS}
setup -B gsl v2_4  -q${build}

# The build system.
setup -B scons v3_0_1  -q p2714b

# The debugger
setup -B gdb v8_0_1

# satellite releases run this setup, then add itself to the following

# where to search for geometry and other configuration
export MU2E_SEARCH_PATH=$MU2E_BASE_RELEASE:$MU2E_DATA_PATH
# Search path for fcl files (overwrites any path from products above)
export FHICL_FILE_PATH=${MU2E_BASE_RELEASE}
# other paths needed to run from this release
export LD_LIBRARY_PATH=`dropit -p $LD_LIBRARY_PATH -sf $MU2E_BASE_RELEASE/lib`
export PYTHONPATH=`dropit -p $PYTHONPATH -sf $MU2E_BASE_RELEASE/scripts/build/python`
export PATH=`dropit -p $PATH -sf $MU2E_BASE_RELEASE/bin`
# root6 needs a path to include files to prevent some runtime warnings
export ROOT_INCLUDE_PATH=`dropit -p $ROOT_INCLUDE_PATH -sf $MU2E_BASE_RELEASE`

# Environment variables used by the test build system.
export PACKAGE_SOURCE=${MU2E_BASE_RELEASE}
export BUILD_BASE=${MU2E_BASE_RELEASE}


