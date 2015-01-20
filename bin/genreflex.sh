#!/bin/bash
#
# Build _dict.cpp and _map.cpp files from the classes.h and classes_def.xml files.
# Both targets are built by a single call to this script.
#
# This script is designed to be called from SConstruct.
#
# Inputs:
# $1 - the relative path to the classes.h file
# $2 - the relative path to either the _dict.cpp or the _map.cpp file that needs to be made
# $3 - the set of -I arguments for the genreflex command
#
# Step 1:
#  - Parse the arguments and build the local shell variables needed to do the work.
#
# Step 2:
#  - Run genreflex.  This creates the _dict.cpp file and a file named classes_ids.cpp.
#
# Step 3:
#  - Rename classes_ids.cc to _map.cpp
#
# Notes:
# 1) The original code built both targets regardless of which of the two was specified
#    As best I can tell scons always tries to build the _dict.cpp file first.  So
#    we could simplify the code to presume that always to be true.  But I cannot guaranttee
#    that behaviour and I left the code with the ability to accept either target as
#    the input argument.
#

usage(){
 echo " "
 echo "Usage: genreflex.sh path_to_classes.h path_to_dict.cpp include_list"
 echo " "
}

if [ "$#" != "3" ]; then
  echo "Illegal arguments to genreflex.sh: there must be exaclty three arguments"
  usage
  exit 2
fi

sourceFile=$1
target=$2
inc=$3

# Check the second argument.  It must be either a relative pathname
# ending in either _dict.cpp or _map.cpp
t1=`expr ${target} : '\(.*\)_dict.cpp'`
t2=`expr ${target} : '\(.*\)_map.cpp'`

if [ -z "$t1" ]; then
 if [ -z "$t2" ]; then
  echo "Illegal argument to the genreflex.sh script. The second argument must be a"
  echo "relative path name ending in either _dict.cpp or _map.cpp"
  echo "The actual value was: " ${target}
  usage
  exit 3
 fi
fi

# Build the names of the two output files, _dict.cpp and _map.cpp
if [ -n "$t1" ];then
  mapFile=${t1}_map.cpp;
  dictFile=${t1}_dict.cpp;
  idsFile=`dirname ${t1}`/classes_ids.cc
elif [ -n "$t2" ]; then
  mapFile=${t2}_map.cpp;
  dictFile=${t2}_dict.cpp;
  idsFile=`dirname ${t2}`/classes_ids.cc
fi;

# Define the flats needed to run the genreflex command.
flags1="--deep --fail_on_warnings --iocomments --capabilities=classes_ids.cc"
flags2=" -D_REENTRANT -DGNU_SOURCE -DGNU_GCC -D__STRICT_ANSI__"
flags3=" -DPROJECT_NAME=\"mu2e\" -DPROJECT_VERSION=\"development\""
flags="${flags1} ${flags2} ${flags3}"

# Build the name of the classes_def.xml input file
classesDef=`dirname ${sourceFile}`/classes_def.xml

# Steps 2 and 3
genreflex ${sourceFile} -s ${classesDef} ${inc} -o $dictFile ${flags}
mv ${idsFile} $mapFile;
