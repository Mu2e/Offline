# Configuration file for Example/03
#
# $Id: ex03.py,v 1.1 2009/09/30 22:57:47 kutschke Exp $
# $Author: kutschke $
# $Date: 2009/09/30 22:57:47 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("Ex03")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Read events from a file (made by example 2)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("ExampleData/ex01Input_1.root")
)

# Look at provenance information.
process.provinspect = mu2e.EDAnalyzer(
    "Ex03InspectProvenance"
    )

# Adjust configuration of message logger.
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.categories.append("ProvenanceInfo")
process.MessageLogger.cerr.Provenance = mu2e.untracked.PSet(
    limit = mu2e.untracked.int32(-1)
)


# Tell the system in which order to execute the modules.
# It is implict that the source module is first.
process.p = mu2e.EndPath( process.provinspect )
