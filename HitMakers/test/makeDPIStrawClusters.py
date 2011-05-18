# Configuration file for making StrawHits.
#
# $Id: makeDPIStrawClusters.py,v 1.5 2011/05/18 02:27:16 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:16 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("HitTest01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("makeDPIStrawClusters.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, they are scheduled later.

# Read events from a file (made by Mu2eG4 example g4test_03.py)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)
process.ConditionsService = mu2e.Service("ConditionsService",
      conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)
# Form StrawHits (SH).
process.makeSH = mu2e.EDProducer(
    "MakeStrawHit",
    g4ModuleLabel = mu2e.string("g4run"),
    seed=mu2e.untracked.vint32(7790),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)
# Check the StrawHits.
process.testSH = mu2e.EDAnalyzer("ReadStrawHit",
    makerModuleLabel = mu2e.string("makeSH"),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# make  the DPIStrawClusters.
process.makeSC = mu2e.EDProducer(
    "MakeDPIStrawCluster",
    makerModuleLabel = mu2e.string("makeSH"),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Check the StrawClusters.
process.testSC = mu2e.EDAnalyzer("ReadDPIStrawCluster",
       g4ModuleLabel = mu2e.string("g4run"),
    makerModuleLabel = mu2e.string("makeSH"),
    clmakerModuleLabel = mu2e.string("makeSC"),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)



# Write an output file.
#process.outfile = mu2e.OutputModule(
#    "PoolOutputModule",
#    fileName = mu2e.untracked.string('file:hits_03.root'),
#    fastCloning = cms.untracked.bool(False),
#)


# End of the section that defines and configures modules.

# Tell the system to execute all paths.
#process.output = mu2e.EndPath(  process.makeSH*process.testSH*process.makeSC*process.testSC*process.outfile );
process.output = mu2e.EndPath(  process.makeSH*process.testSH*process.makeSC*process.testSC );

