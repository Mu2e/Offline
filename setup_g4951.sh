#
# $Id: setup_g4951.sh,v 1.1 2012/04/19 22:58:09 genser Exp $
# $Author: genser $
# $Date: 2012/04/19 22:58:09 $
#
# Original author KLG based on setup.sh
#
# Interim setup to introduce geant4 9.5.p01

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

getfn() { fn="$1"; while [ -L "$fn" ]; do fn="`readlink $fn`"; done; echo "$fn"; }

# This variable contains the absolute path to the directory
# that contains this file  (regardless of cwd when this script is sourced ).
myLocation=`cd "$(dirname $(getfn ${BASH_SOURCE}) )" >/dev/null 2>&1 && echo \$PWD`

source ${myLocation}/setup.sh

# Setup the framework and its dependent products
setup art v1_00_11 -qmu2e:prof
setup splines v1_00_01 -q a7:prof

# Geant4 and its cross-section files.
setup geant4 v4_9_5_p01 -qgcc46:prof

setup g4abla v3_0
setup g4emlow v6_23
setup g4neutron v4_0
setup g4neutronxs v1_1
setup g4photon v2_2
setup g4pii v1_3
setup g4radiative v3_4
setup g4surface v1_0

# A hack that we hope can go away soon.
export G4LIBDIR=$G4LIB/$G4SYSTEM
G4INCLUDE=${G4INCLUDE}/Geant4
