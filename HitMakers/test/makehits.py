
# Configuration file for Readback
#
# $Id: makehits.py,v 1.10 2010/08/18 23:14:03 logash Exp $
# $Author: logash $
# $Date: 2010/08/18 23:14:03 $
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
    input = mu2e.untracked.int32(200)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("makehits.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Read events from a file (made by Mu2eG4 example g4test_03.py)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Form CrudeStrawHits (CSH).
process.makeCSH = mu2e.EDProducer(
    "MakeStrawHit",
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Check the crudeStrawHits.
process.testCSH = mu2e.EDAnalyzer("ReadStrawHit",
    diagLevel    = mu2e.untracked.int32(3),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Make some clusters.
process.clustertest = mu2e.EDAnalyzer("ClusterHackv00",
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Write an output file.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:hits_03.root'),
    fastCloning = cms.untracked.bool(False),
)


# End of the section that defines and configures modules.

# Tell the system to execute all paths.
#process.output = mu2e.EndPath(  process.makeCSH*process.testCSH*process.clustertest*process.outfile );
process.output = mu2e.EndPath(  process.makeCSH*process.testCSH*process.outfile );


