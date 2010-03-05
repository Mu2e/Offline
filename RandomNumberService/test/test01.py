# Configuration file to test RandomNumberGenerator
#
# $Id: test01.py,v 1.1 2010/03/05 16:07:38 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/03/05 16:07:38 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("RNGTest01")

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
                          globalSeed=mu2e.untracked.int32(9877)
))

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

process.source  = mu2e.Source("EmptySource")
process.rngtest = mu2e.EDAnalyzer("RNGTest")
#process.rngtest = mu2e.EDAnalyzer("RandomNumberSaver")

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.rngtest );

