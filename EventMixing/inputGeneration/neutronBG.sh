#! /bin/bash
#
# A grid script to generate neutron backgrounds in the tracker and calorimeter.
#
#   $Id: neutronBG.sh,v 1.1 2012/07/15 16:01:19 kutschke Exp $
#   $Author: kutschke $
#   $Date: 2012/07/15 16:01:19 $
#
# It takes 4 arguments:
#   cluster    - the condor job cluster number
#   process    - the condor job process number, within the cluster
#   user       - the username of the person who submitted the job
#   submitdir  - the directory from which the job was submitted ( not used at this time).
#
# Outputs:
#  - All output files are created in the grid working space.  At the end of the job
#    all files in this directory will be copied to:
#      /grid/data/mu2e/outstage/$user/${cluster}_${process}
#    This includes a copy of any input files that the user does not first delete.
#
# Notes:
#
# 1) For documentation on using the grid, see
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml
#    For details on cpn and outstage see:
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#cpn
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#outstage
#
# 2) To test if there are files to move from the original current working
#    directory to the self destructing working directory, we use
#
#      if (( $( ls -1 $ORIGDIR | wc -l) > 0 )); then ... ; fi
#
#    This is potentially fragile.  If a user aliases or otherwise redefines
#    ls, this may fail.  Unfortunately the grid worker nodes do not let
#    us write /bin/ls instead of plain ls.  If we do, then we get an
#    error
#
#       ls: write error: Broken pipe
#
#    So we are living with the fragile code.
#

# Copy arguments into meaningful names.
cluster=$1
process=$2
user=$3
submitdir=$4
outstagebase=$5
sprocess=`printf "%4.4d" $process`
echo "Input arguments:"
echo "Cluster:    " $cluster
echo "Process:    " $process
echo "User:       " $user
echo "SubmitDir:  " $submitdir
echo "Outstage base directory: " $outstagebase
echo "Padded Process:          " $sprocess
echo " "
echo "Hostname:   " `hostname`
echo " "

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
ORIGDIR=`pwd`
TMP=`mktemp -d ${OSG_WN_TMP:-/var/tmp}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${OSG_WN_TMP:-/var/tmp}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0
cd $TMP

# If there are files in the original grid working directory, move them to the new directory.
# Fragile!!! See note 2.
if (( $( ls -1 $ORIGDIR | wc -l) > 0 )); then
 mv $ORIGDIR/* .
fi
# End of the section you should not change.

# Directory in which to put the output files.
outstage=${outstagebase}/$user

# Establish environment.
source /grid/fermiapp/products/mu2e/setupmu2e-art.sh
source /grid/fermiapp/mu2e/users/kutschke/HEAD_120711/Offline_g4951/setup_g4951.sh
source /grid/fermiapp/mu2e/users/kutschke/HEAD_120711/Offline_g4951/bin/addlocal.sh

# Copy the stopping muon position file:
/grid/fermiapp/minos/scripts/cpn $MU2E_DATA_PATH/ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt .

# Set the random number seeds
baseRunNumber=1
maxUniqueEngines=20

# Update the .fcl file to make it specific to this process within this cluster.
updateFCLFile.sh $baseRunNumber $maxUniqueEngines $cluster $process neutronBG.fcl thisjob.fcl

# Run the Offline job.
time mu2e -c thisjob.fcl -n 7000000 >& thisjob.log

# Remove any files that should not be copied to the output staging area.
# Subdirectories will not be copied so there is no need to delete them.
rm neutronBG.fcl
rm stoppedMuons_02.txt

# Create a directory in the output staging area.
/grid/fermiapp/mu2e/bin/createOutStage.sh ${cluster} ${sprocess} ${user} ${outstagebase}

# Copy all files from the working directory to the new directory in the output staging area.
# This will not copy subdirectories.
/grid/fermiapp/minos/scripts/cpn * ${outstage}/${cluster}/${cluster}_${sprocess}

exit 0
