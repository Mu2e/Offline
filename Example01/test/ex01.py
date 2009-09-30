# Configuration file for Example/01
#
# $Id: ex01.py,v 1.1 2009/09/30 22:57:47 kutschke Exp $
# $Author: kutschke $
# $Date: 2009/09/30 22:57:47 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("Ex01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(300)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("ex01histo.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Read events from a file (made by example 2)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("ExampleData/ex01Input_1.root")
)

# Get the hits out of the event, make plots and printout.
process.hitinspect = mu2e.EDAnalyzer(
    "Ex01InspectHits",
    maxFullPrint=mu2e.untracked.int32(10)
    )

# End of the section that defines and configures modules.

# Adjust configuration of message logger.
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.categories.append("ToyHitInfo")
process.MessageLogger.cerr.ToyHitInfo = mu2e.untracked.PSet(
    limit = mu2e.untracked.int32(-1)
)

# Tell the message logger summary to include INFO messages.
process.MessageLogger.cerr_stats.threshold = mu2e.untracked.string('INFO')

# Tell the system in which order to execute the modules.
# It is implict that the source module is first.
process.p = mu2e.EndPath( process.hitinspect )
