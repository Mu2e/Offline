# Configuration file for G4Test03
#
# $Id: houghtest.py,v 1.7 2011/05/18 02:27:16 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:16 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("HoughTest03")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(100)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("houghtest_03.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=cms.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mu2e.Source("EmptySource")

#  Make some generated tracks and add them to the event.
#1 is nothing but conversions
#2 has cosmics and radiative pi's

process.generate = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/genconfigReflection.txt")
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    seed=mu2e.untracked.vint32(9877),
)

# Save state of random numbers to the event.
process.randomsaver = mu2e.EDAnalyzer("RandomNumberSaver")

process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('file:hough_03.root'),
        outputCommands = cms.untracked.vstring(
     'keep *_*_*_*',
#     'drop mu2eSimParticles_*_*_*'   # Uncomment this line to reduce file size.
     ),
)

# Look at the hits from G4.
process.checkhits = mu2e.EDAnalyzer(
# this line tells me where to get the plugin file from
    "HoughTest",
    NPeakSearch = mu2e.uint32(50),
    UseMCHits = mu2e.bool(True),
    hitCreatorName = mu2e.string("g4run"),
    maxFullPrint = mu2e.untracked.int32(5)
)

# End of the section that defines and configures modules.

# Adjust configuration of message logger.
# Enable debug printout from the module instance "hitinspect".
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.cerr.threshold = mu2e.untracked.string('DEBUG')
process.MessageLogger.debugModules.append("hitinspect")
process.MessageLogger.categories.append("ToyHitInfo")
process.MessageLogger.categories.append("GEOM")

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.generate*process.g4run*process.checkhits*process.outfile );

