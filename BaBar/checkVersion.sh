#! /bin/bash
#
# Check the version of the checked out BaBar svn code and compare
# it to the required version, which is passed in as an argument.
#
# The check is done on a single directory.  I believe that this is enough
# but I can't be certain that it won't change in the future.
#
# This code only prints warnings; it does not do an update.
#
# $Id: checkVersion.sh,v 1.1 2013/04/10 16:01:25 kutschke Exp $
# $Date: 2013/04/10 16:01:25 $
# $Author: kutschke $
#
# Contact person Rob Kutschke
#

# The required version number is given as an input argument.
requiredVersion=$1

if ! [[ -f "SConstruct" ]]; then
  echo "ERROR: must be executed from top level of offline code." 1>&2
  exit 1
fi

# Holds the output of an svn status command.
tmpFile="svntmp"

if [ -e $tmpFile ]; then
  /bin/rm -f $tmpFile
  echo  "remove it"
fi

cmd="svn status -u -q BaBar/BaBar"
$cmd >& $tmpFile

# Check that the output of the svn command looks correct.
nLines=`wc $tmpFile | awk '{print $1}'`

if [ $nLines != 1 ]; then
  echo "ERROR: Expected exactly one line from the command " $cmd
  echo "Number of lines is: " $nLines
  echo "File contents: "
  echo " "
  cat $tmpFile
  echo " "
  echo "Please remove the file " $tmpFile " when you are done checking this."
  exit 1
fi

nLines2=`cat $tmpFile | grep "Status against revision:" | wc | awk '{print $1}'`
if [ $nLines != 1 ]; then
  echo "ERROR: Did not find the expected line in the output of the command " $cmd
  echo "File contents: "
  echo " "
  cat $tmpFile
  echo " "
  echo "Please remove the file " $tmpFile " when you are done checking this."
  exit 1
fi

checkedOutVersion=`cat $tmpFile | grep "Status against revision:" | awk '{print $4}' `

if [ $checkedOutVersion != $requiredVersion ]; then
  echo "The checked out version of the BaBar code does not match the required version."
  echo "Checked out version: " $checkedOutVersion
  echo "Required version:    " $requiredVersion
  echo "If you need to update the BaBaR code you should: source BaBar/update.sh"
fi

if [ $checkedOutVersion == $requiredVersion ]; then
  echo "The checked out BaBar code matches the required version: " $requiredVersion
fi

/bin/rm -f $tmpFile
