# Drop all information except the saved random number engine state.
#
# $Id: drop.py,v 1.3 2010/09/27 20:01:46 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/09/27 20:01:46 $
#
# Original author Rob Kutschke
#
# Notes:
# 1) There is a bug in the interaction between root and the persistency
#    mechanism in the framework.  The problem is that fast cloning does
#    not work.  The solution is to disable fast cloning, which is done
#    in the configuration of the output modules.
#

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("Drop01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")


# Start each new event with an empty event.
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root"),
)

# Define the output file. See note 1.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:drop_03.root'),
    fastCloning = mu2e.untracked.bool(False),
    outputCommands = cms.untracked.vstring(
     'drop *_*_*_*',
     'keep edmRNGsnapshots_*_*_*'
     )
)

# End of the section that defines and configures modules.


# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.outfile );
