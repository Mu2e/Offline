#!/bin/bash
#
# Build _dict.cpp file from the classes.h and classes_def.xml files.
#
# This script is designed to be called from SConstruct.
#
# Inputs:
# $1 - Source: the relative path to the classes.h file
# $2 - Target: the relative path to either the _dict.cpp file that needs to be made
# $3 - the set of -I arguments for the genreflex command
# $4 - the name of the _dict.so file that will be built from 2) in a later step
#
# Step 1:
#  - Parse the arguments and build the local shell variables needed to do the work.
#
# Step 2:
#  - Run genreflex.  This creates the _dict.cpp file in the temp directory.
#    It also creates 2 files in the lib directory: the _rdict.pcm and .rootmap files.
#

usage(){
 echo " "
 echo "Usage: genreflex.sh path_to_classes.h path_to_dict.cpp \\"
 echo "                     include_list lib_name prof_or_debug"
 echo " "
}

if [ "$#" != "5" ]; then
  echo "Illegal arguments to genreflex.sh: there must be exactly five arguments"
  usage
  exit 2
fi

# Process the argument list.
sourceFile=$1
dictFile=$2
inc=$3
libname=$4
debug_level=$5

# Define the flags needed to run the genreflex command.
flags="--fail_on_warnings"
if [ "${debug_level}" == "prof" ]; then
 flags=${flags}" -DNDEBUG"
fi

# Build the name of the classes_def.xml input file
# A previous step checked for its existence.
classesDef=`dirname ${sourceFile}`/classes_def.xml

# Build the name of the rootmap file
rootmapFile=`echo $libname | sed 's/_dict\.so/_dict\.rootmap/'`

# Step 2
genreflex ${sourceFile} -s ${classesDef} ${inc} -l $libname -o $dictFile ${flags} \
--rootmap-lib=$libname \
--rootmap=$rootmapFile

unset sourceFile
unset dictFile
unset inc
unset libname
unset debug_level
unset classesDef
unset rootmapFile
