#! /bin/bash
#
#
# Run this script from within a grid job to append per job information to
# the end of the .fcl file:
#   - the run number for this job
#   - the baseSeed for the random number engines in this job.
#   - the maximum number of random number engines that may be seeded in  this job.
#
# See also:
#   http://mu2e.fnal.gov/atwork/computing/RandomBasics.shtml
#   http://mu2e.fnal.gov/atwork/computing/Random.shtml
#
# Arguments:  all are input only and all are required.
#   1 - the base run number for this grid cluster
#   2 - the maximum number of unique random engines per job
#   3 - the grid cluster number ( reserved for future use; not used in this script ).
#   4 - the grid process number
#   5 - the name of the input fcl file
#   6 - the name of the output fcl file; it contains the
#
# What this script does:
#
# 1) It copies the input fcl file to the output file.
#
# 2) It appends a line to the output fcl file to set the run number to the base run number
#    plus the process number.  By art convention, run numbers start at 1, not 0.
#
# 3) This script collaborates with the SeedService to implement the following algorithm.
#    For each job the base seed is set by this script to
#       baseSeed=$(( $runNumber * $maxUniqueEngines ))
#    The values of baseSeed and maxUniqueEngines are passed to the SeedService by appending
#    their values to the output file.
#
#    The SeedService will generate unique seeds in the half open range:
#       [baseSeed, baseSeed+maxUniqueEngines )
#
#    If this does not provide enough seeds for the job, the SeedService will throw an
#    exception and stop the job.  The exception will occur during the initialization phase
#    of the job, in either the c'tor or beginRun phase.
#
#    To recover from this error, you should increase the value of maxUniqueEngines and resubmit your
#    job.
#
# One can image that future versions of this script that might cooperate with the SeedService in
# some other way.  Therefore the concept of the base seed is distinguished from the run number,
# even though, in this implementation, they are tightly coupled.
#

if [ $# != 6 ]; then
  echo "updateFCLFile.sh: You must supply a exactly six arguments."
  echo "Number of arguments: " $#
  exit 1
fi

# Copy the arguments into useful names.
baseRunNumber=$1
maxUniqueEngines=$2
cluster=$3
process=$4
inputFile=$5
outputFile=$6

# Safety checks
if [ $baseRunNumber -lt 1 ]; then
  echo "updateFCLFILE: baseRunNumber must be greater than 0;  it is: " $baseRunNumber
  exit 1
fi

if [ ! -e $inputFile ]; then
  echo "updateFCLFILE: the input file," $inputFile "does not exist."
  exit 1
fi

if [ $maxUniqueEngines -lt 1 ]; then
  echo "updateFCLFILE: the maxUniqueEngines parameter must be > 1; it is" $maxUniqueEngines
  exit 2
fi

# Compute the derived quantities.
runNumber=$(( $baseRunNumber + $process ))
baseSeed=$(( ( $runNumber - 1 ) * $maxUniqueEngines ))


# Write the output file:
cp $inputFile $outputFile
echo "source.firstRun                            : " $runNumber        >> $outputFile
echo "services.user.SeedService.baseSeed         : " $baseSeed         >> $outputFile
echo "services.user.SeedService.maxUniqueEngines : " $maxUniqueEngines >> $outputFile
