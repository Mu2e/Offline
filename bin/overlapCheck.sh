#! /bin/bash
#
# Run the root-based check for overlapping geometry objects in a gdml file
#  - also checks for errors in the loading of the gdml
#

usage() {
cat <<EOF

   overlapCheck.sh [OPTIONS] [GDMLFILE]

   run the root overlap check, return error code if the process failed
   or if overlaps were detected.
   if GDMLFILE argument is not proesent, look in the Muse standard location
   \$MUSE_BUILD_DIR/Offline/gen/gdml/muse.gdml.

   -q don't print the full output, only errors

   -b if no GDMLFILE provided, try to build it in Muse standard location
      then run checks on it

EOF
}

# the number of lines normally printed by the root command
NEXPECTED=17

QUIET=
BUILD=
REPORT=

while getopts hqbr OPTION ; do
    case $OPTION in
        h) usage; exit 0 ;;
        q) QUIET=yes ;;
        b) BUILD=yes ;;
        r) REPORT=yes ;;
        *) usage; exit 1 ;;
    esac
done
shift $(($OPTIND-1))
GDMLFILE="$@"



if [ -z "$GDMLFILE" ]; then
    GDMLFILE=$MUSE_BUILD_DIR/Offline/gen/gdml/mu2e_common.gdml
    if [[ ! -e "$GDMLFILE" && -n "$BUILD" ]]; then
        muse build GDML
    fi
fi

if [ ! -e "$GDMLFILE" ]; then
    echo "ERROR - no gdml file specified or unable to find/build it"
    exit 1
fi


SCRIPTDIR=$( dirname $(readlink -f $BASH_SOURCE) )
LL=$( root -l -b -q $SCRIPTDIR/overlapCheck.C\(\"$GDMLFILE\"\) 2>&1  )
RRC=$?

[ -z "$QUIET" ] && echo "$LL"

NLINES=$( echo "$LL" | wc -l )
NERROR=$( echo "$LL" | grep -ci error )

NILLEGAL=$( echo "$LL" | grep "illegal" | awk '{print $NF}' )
[ -z "$NILLEGAL" ] && NILLEGAL=0

if [ -n "$REPORT" ]; then
    TOTALS=$( echo "$LL" | grep "in Geometry imported from GDML" | \
        awk '{for(i=4;i<9;i++) printf "%s ", $i }')
    echo "Overlaps result: $TOTALS, $NILLEGAL illegal"
    if [ $NILLEGAL -gt 0 ]; then
       echo "$LL" | awk 'BEGIN{flag=0;}{if(flag==1) print $0; if($1=="===") flag=1; }'
    fi
fi

RC=0
if [ $RRC -ne 0 ]; then
    echo "ERROR - root command failed"
    RC=1
fi
if [ $NERROR -ne 0 ]; then
    echo "ERROR - error messages from root"
    RC=1
fi
if [[ $NLINES -ne $NEXPECTED && $NILLEGAL -eq 0 ]]; then
    echo "ERROR - unexpected number of lines in root response: $NLINES"
    RC=1
fi
if [ $NILLEGAL -ne 0 ]; then
    echo "ERROR - overlaps found"
    RC=1
fi

exit $RC
