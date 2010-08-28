#
# $Id: setup.sh,v 1.12 2010/08/28 18:31:16 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/08/28 18:31:16 $
#
# Original author Rob Kutschke
#
# Setup the environment to build a full release of the Mu2e software.
# This presumes that you have already established the Mu2e environment
# for this machine with something like:
# source /prj/mu2e/

if [ "`basename $0 2>/dev/null`" = "setup.sh" ];then
    echo "You should be sourcing this file, not executing it."
    exit
fi

if [ "${MU2E}" = '' ];then
    echo "The environment variable MU2E is not set."
    echo "You must setup the local Mu2e environment before sourcing this script."
    exit
fi

# Define the directory in which this file lives as the root of a release.
export MU2E_BASE_RELEASE=`cd "$(dirname ${BASH_SOURCE})" >/dev/null 2>&1 && echo \$PWD`
echo "Base release directory is: " $MU2E_BASE_RELEASE

# This will setup all products on which framework depends.
setup framework v1_1_3

# Chose version of G4 and its cross-section files.
setup geant4 v4_9_3_p01 -q g77-OpenGL
setup g4neutron   v3_13a
setup g4emlow     v6_2
setup g4photon    v2_0
setup g4radiative v3_2
setup g4abla      v3_0

# Other products
setup heppdt      v3_04_01

# Tell the framework to look in the local area to find modules.
source ${MU2E_BASE_RELEASE}/bin/setup_mu2e_project.sh

# Build the symlink directories for the BaBar code.
# Only do so if the BaBar package is checked out locally.
if [  -f "./BaBar/makeInclude.csh" ]
then
  ./BaBar/makeInclude.csh
fi
