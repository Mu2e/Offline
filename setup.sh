# Setup the environment to build or use a full release of the Mu2e software.
# This checks that you have already established the Mu2e environment.
#
# Original author Rob Kutschke
#
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

#
# if this is a partial checkout, run the base release setup
# and add this to the path
#

# This variable contains the physical path to the directory
# that contains this file  (regardless of cwd when this script is sourced ).
SCRIPTDIR=`dirname $(readlink -f $BASH_SOURCE)`
BASEREPO=$( git config --file $SCRIPTDIR/.git/config --get mu2e.baserelease )

if [ -n "$BASEREPO" ]; then

  #echo "found partial checkout active"
  #echo "base release: $BASEREPO"
  if [ ! -f $BASEREPO/setup.sh ]; then
    echo " ERROR - could not find base release setup"
    return 1
  fi

  # do the base release setup (redefines SCRIPTDIR)
  source $BASEREPO/setup.sh

  # this file
  SATSCRIPTDIR=`dirname $(readlink -f $BASH_SOURCE)`
  # Add the satellite release to path variables.
  export MU2E_SATELLITE_RELEASE=$SATSCRIPTDIR
  export MU2E_SEARCH_PATH=`dropit -p $MU2E_SEARCH_PATH -sf $MU2E_SATELLITE_RELEASE`
  export FHICL_FILE_PATH=`dropit -p $FHICL_FILE_PATH -sf $MU2E_SATELLITE_RELEASE`
  export CET_PLUGIN_PATH=`dropit -p $CET_PLUGIN_PATH -sf $MU2E_SATELLITE_RELEASE/lib`
  export LD_LIBRARY_PATH=`dropit -p $LD_LIBRARY_PATH -sf $MU2E_SATELLITE_RELEASE/lib`
  export PYTHONPATH=`dropit -p $PYTHONPATH -sf $MU2E_SATELLITE_RELEASE/scripts/build/python`
  export PATH=`dropit -p $PATH -sf $MU2E_SATELLITE_RELEASE/bin`
  export ROOT_INCLUDE_PATH=`dropit -p $ROOT_INCLUDE_PATH -sf $MU2E_SATELLITE_RELEASE`

  if [ -f $MU2E_SATELLITE_RELEASE/.buildopts ] ; then
      echo "WARNING - buildopts will be taken from the base release, "
      echo "        in partial checkout your local buildopts is ignored"
  fi
  # items below were defined by the base release
  # source setup.sh, so exit now
  return 0

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
export MU2E_UPS_QUALIFIERS=+e20:+${build}
export MU2E_ART_SQUALIFIER=s108

MU2E_G4_GRAPHICS_QUALIFIER=''
if [[ $($MU2E_BASE_RELEASE/buildopts --g4vis) == qt ]]; then
    MU2E_G4_GRAPHICS_QUALIFIER=':+qt'
fi

MU2E_G4_VECGEOM_QUALIFIER=''
if [[ $($MU2E_BASE_RELEASE/buildopts --g4vg) == on ]]; then
    MU2E_G4_VECGEOM_QUALIFIER=':+vg'
fi

MU2E_G4_MT_QUALIFIER=''
if [[ $($MU2E_BASE_RELEASE/buildopts --g4mt) == on ]]; then
    MU2E_G4_MT_QUALIFIER=':+mt'
fi

export MU2E_G4_EXTRA_QUALIFIER=''

# Setup the framework and its dependent products
setup -B art v3_09_00 -q${MU2E_UPS_QUALIFIERS}
setup -B art_root_io v1_08_00 -q${MU2E_UPS_QUALIFIERS}

# Geant4 and its cross-section files.
if [[ $($MU2E_BASE_RELEASE/buildopts --trigger) == "off" ]]; then
  setup -B geant4 v4_10_7_p01c -q${MU2E_UPS_QUALIFIERS}${MU2E_G4_GRAPHICS_QUALIFIER}${MU2E_G4_VECGEOM_QUALIFIER}${MU2E_G4_MT_QUALIFIER}${MU2E_G4_EXTRA_QUALIFIER}
else
  setup -B xerces_c v3_2_3   -q${MU2E_UPS_QUALIFIERS}
fi

# Get access to raw data formats.
setup -B mu2e_artdaq_core v1_05_09_01 -q${MU2E_UPS_QUALIFIERS}:+${MU2E_ART_SQUALIFIER}

setup -B heppdt   v03_04_02 -q${MU2E_UPS_QUALIFIERS}
setup cetpkgsupport
setup -B KinKal   v00_01_06a  -q${MU2E_UPS_QUALIFIERS}:p392
setup -B BTrk   v1_02_31  -q${MU2E_UPS_QUALIFIERS}:p392
setup -B cry   v1_7n  -q${MU2E_UPS_QUALIFIERS}
setup -B gsl v2_6a
setup curl v7_64_1
setup cryptopp v08_02_00 -q${MU2E_UPS_QUALIFIERS}

# The build system.
setup -B scons v3_1_2a  -q +p392

# The debugger
setup -B gdb v9_2

# satellite releases run this setup, then add itself to the following

# where to search for geometry and other configuration
export MU2E_SEARCH_PATH=$MU2E_BASE_RELEASE:$MU2E_DATA_PATH
# Search path for fcl files (overwrites any path from products above)
export FHICL_FILE_PATH=${MU2E_BASE_RELEASE}
# other paths needed to run from this release
export CET_PLUGIN_PATH=`dropit -p $CET_PLUGIN_PATH -sf $MU2E_BASE_RELEASE/lib`
export LD_LIBRARY_PATH=`dropit -p $LD_LIBRARY_PATH -sf $MU2E_BASE_RELEASE/lib`
export PYTHONPATH=`dropit -p $PYTHONPATH -sf $MU2E_BASE_RELEASE/scripts/build/python`
export PATH=`dropit -p $PATH -sf $MU2E_BASE_RELEASE/bin`
# root6 needs a path to include files to prevent some runtime warnings
export ROOT_INCLUDE_PATH=`dropit -p $ROOT_INCLUDE_PATH -sf $MU2E_BASE_RELEASE`

