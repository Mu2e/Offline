
# Configuration file for reading back StrawHits from a file.
#
# $Id: read.py,v 1.6 2010/08/26 19:16:45 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/08/26 19:16:45 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("ReadTest01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(200)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("readhits.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Read events from a file (made by Mu2eG4 example g4test_03.py)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("hits_03.root")
)

# Check the crudeStrawHits.
process.testSH = mu2e.EDAnalyzer("ReadStrawHit",
    makerModuleLabel = mu2e.string("makeSH"),
    diagLevel    = mu2e.untracked.int32(3),
    maxFullPrint = mu2e.untracked.int32(100)
)

# End of the section that defines and configures modules.

# Tell the system to execute all paths.
process.output = mu2e.EndPath( process.testSH );

