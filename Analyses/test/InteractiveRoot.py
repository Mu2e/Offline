# Configuration file for testing interactive root.
#
# $Id: InteractiveRoot.py,v 1.2 2010/09/24 17:02:29 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/09/24 17:02:29 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("InteractiveRoot")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(100)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("InteractiveRoot.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Read events from a file (made by example 3)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Look at the hits from G4.
process.root1 = mu2e.EDAnalyzer(
    "InteractiveRoot",
    g4ModuleLabel = mu2e.string("g4run"),
)

# End of the section that defines and configures modules.

# Tell the system to execute the modules.
process.output = mu2e.EndPath(  process.root1 );

