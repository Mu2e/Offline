
# This configuration file is sourced from 
# https://github.com/Mu2e/codetools/blob/master/bin/github/jenkins_tests/mu2e-offline-build-test/build.sh

# These environment variables influence which art jobs which must return status 0 for a build test to pass.
# For each array index, the elements from each of the 3 arrays 
# are combined as below:
#  mu2e -n "${NEVTS_TJ[$i]}" -c "${FCLFILES[$i]}" >> "${JOBNAMES[$i]}.log" 2>&1
# where $i is the index corresponding to an element of each of the 3 arrays.

# an easily identifiable, human-recognisable label for the art job being run
declare -a JOBNAMES=( \
  "ceSimReco" \
  "g4test_03MT" \
  "transportOnly" \
  "PS" \
  "g4study" \
  "cosmicSimReco")

# a path to the fcl file that is run in the test
declare -a FCLFILES=("Validation/fcl/ceSimReco.fcl" \
  "Mu2eG4/fcl/g4test_03MT.fcl" \
  "Mu2eG4/fcl/transportOnly.fcl" \
  "JobConfig/beam/PS.fcl" \
  "Mu2eG4/fcl/g4study.fcl" \
  "Validation/fcl/cosmicSimReco.fcl")

# the number passed to the art -n flag
declare -a NEVTS_TJ=("1" "10" "1" "1" "1" "1")
