#
# $Id: hitDisplay.py,v 1.1 2011/05/17 00:28:23 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/17 00:28:23 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("HitDisplay")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(5)
)

# Define the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Define the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("hitDisplay.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define the random number generator service.
#process.RandomNumberGeneratorService = mu2e.Service("RandomNumberGeneratorService")

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# Read events from a file (made by example 3)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Look at the hits from G4.
process.hitDisplay = mu2e.EDAnalyzer(
    "HitDisplay",
    generatorModuleLabel = mu2e.string("generate"),
    g4ModuleLabel        = mu2e.string("g4run"),
    hitMakerModuleLabel  = mu2e.string("makeSH"),
    trackerStepPoints    = mu2e.string("tracker"),
    minEnergyDep         = mu2e.double(0.0001),
    minHits              = mu2e.uint32(5),
    doDisplay            = mu2e.untracked.bool(True),
)

# End of the section that defines and configures modules.

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.hitDisplay );
