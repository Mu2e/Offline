# Read the output file from g4test_03.py and rewrite it,
# keeping only the random number state.
#
# $Id: drop03.py,v 1.1 2010/10/13 23:41:51 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/10/13 23:41:51 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("Drop03")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Define the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Start each new event with an empty event.
# Read events from a file (made by example 3)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Define the output file and its minimal contents.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:drop_03.root'),
    outputCommands = mu2e.untracked.vstring(
     'drop *_*_*_*',
     'keep *_randomsaver_*_*'
    ),
)

# End of the section that defines and configures modules.

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.outfile);
