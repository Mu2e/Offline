# Configuration file to test RandomNumberGenerator
#
# $Id: test03.py,v 1.1 2010/03/07 22:01:00 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/03/07 22:01:00 $
#
# Restore state at start of each event
# Test fire instances of Randflat and RandGaussQ on odd
# numbered events; compare to output of test01.py
#
# Original author Rob Kutschke

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("RNGTest03")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(5)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberService",
             restoreStateLabel=mu2e.untracked.string("randomsaver"),
             debug=mu2e.untracked.bool(True)
))

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("randomtest_01.root")
)

process.rngtest     = mu2e.EDAnalyzer("RNGTest",
                       doSkip= mu2e.untracked.int32(2),
                      )
process.outfile     = mu2e.OutputModule("PoolOutputModule",
                      fileName = mu2e.untracked.string('file:randomtest_03.root'),
)

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.rngtest*process.outfile );

