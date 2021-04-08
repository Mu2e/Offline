#! /bin/bash
SCRIPTDIR=$( dirname $(readlink -f $BASH_SOURCE) )
root -l -b << EOF
.x $SCRIPTDIR/overlapCheck.C("$1");
.q
EOF
