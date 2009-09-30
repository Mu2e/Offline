#
# $Id: setups.sh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
# $Author: kutschke $
# $Date: 2009/09/30 22:57:47 $
#
# Original author Rob Kutschke
#
# Make sure that ups is set up.
# This really belongs in your .profile, .shrc or the like.
#

if [  -f "/fnal/ups/etc/setups.sh" ]
then
    . "/fnal/ups/etc/setups.sh"
fi
