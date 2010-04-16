#
# $Id: setup.sh,v 1.4 2010/04/16 14:49:57 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/04/16 14:49:57 $
#
# Original author Rob Kutschke
#
# The stuff in this first file should go into .profile.
# Then remove this line:
source setups.sh

# On ilcsim*

source  /prj/mu2e/setup.sh


# On fnalu, 32 bit sl4.
#source /afs/fnal.gov/files/code/mu2e/d3/mu2e_v0.SL4.i686/bin/setup_mu2e.sh


# Working area
source  bin/setup_mu2e_project.sh

# Geant is not yet in the externals package; get it from ups.
setup geant4 v4_9_2_p01 -q g77-OpenGL

# various G4 cross-section data files

setup g4neutron   v3_13a
setup g4emlow     v6_2
setup g4photon    v2_0
setup g4radiative v3_2
setup g4abla      v3_0

# Other products
setup heppdt      v3_04_01
