#
# $Id: setup.sh,v 1.18 2011/05/17 18:43:57 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/17 18:43:57 $
#
# Original author Rob Kutschke
#
# Setup the environment to build a full release of the Mu2e software.
# This presumes that you have already established the Mu2e environment
# for this machine with something like:
# source /prj/mu2e/

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

case ${EXTERNALSVERSION} in
  2)
    # This will set up G4 also.
    setup art v0_06_03 -qa2:debug

    export FW_HOME=${ART_FQ_DIR}
    export FRAMEWORK_DIR=${FW_HOME}

    setup geant4 v4_9_4_p01 -qgcc45

    setup g4neutron v3_14
    setup g4emlow v6_19
    setup g4photon v2_1
    setup g4radiative v3_3
    setup g4abla v3_0

    setup heppdt v3_04_01 -qgcc45

    setup scons v1_3_0a -qgcc45

    # This is not as expected in the product
    export CLHEP_LIB=${CLHEP_DIR}/lib

    export FHICL_FILE_PATH=.:fcl;

  ;;
  1)
    setup framework v1_1_4
    export FRAMEWORK_DIR=$FW_HOME

    # For historical reasons, Geant4 has different qualifiers on SLF4 and SLF5.
    mu2eG4Qual="gfortran-OpenGL-GDML"
    if grep 'release 4' /etc/redhat-release >/dev/null;then
      mu2eG4Qual=g77-OpenGL-GDML  
    fi

    # Chose version of G4 and its cross-section files.
    setup geant4 v4_9_3_p02 -q $mu2eG4Qual
    
  ;;
  *)
    # Deprecated: backwards compatibility with the old externals system.
    setup framework v1_1_3
    setup geant4 v4_9_3_p01 -q g77-OpenGL        
    setup heppdt      v3_04_01
esac

if (( ${EXTERNALSVERSION:-0} < 2 )); then
  # G4 cross-section files.
  setup g4neutron   v3_13a
  setup g4emlow     v6_2
  setup g4photon    v2_0
  setup g4radiative v3_2
  setup g4abla      v3_0

  # Other products
  setup heppdt      v3_04_01
fi

# Another depracated section: only needed for old externals system.
if [ "${EXTERNALSVERSION}" = '0' ]; then
    export CLHEP_INC=${CLHEP_DIR}/include
    export CLHEP_LIB=${CLHEP_DIR}/lib
    export BOOST_LIB=${BOOST_DIR}/lib
    export HEPPDT_INC=${HEPPDT_DIR}/include
    export HEPPDT_LIB=${HEPPDT_DIR}/lib
    export LIBSIGCPP_INC=${LIBSIGCPP_DIR}/include
    export LIBSIGCPP_LIB=${LIBSIGCPP_DIR}/lib
    export G4LIB=${G4LIB}/Linux-g++
fi

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

# Build the symlink directories for the BaBar code.
# Only do so if the BaBar package is checked out locally.
if [  -f "./BaBar/makeInclude.csh" ]
then
  ./BaBar/makeInclude.csh
fi
