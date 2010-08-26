#
# $Id: setup.sh,v 1.8 2010/08/26 15:50:24 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/08/26 15:50:24 $
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

# This will setup all products on which framework depends.
setup framework v1_1_2

# Chose version of G4 and its cross-section files.
setup geant4 v4_9_3_p01 -q g77-OpenGL
setup g4neutron   v3_13a
setup g4emlow     v6_2
setup g4photon    v2_0
setup g4radiative v3_2
setup g4abla      v3_0

# Other products
setup heppdt      v3_04_01

# Working area
source bin/setup_mu2e_project.sh

# Build the symlink directories for the BaBar code.
./BaBar/makeInclude.csh

