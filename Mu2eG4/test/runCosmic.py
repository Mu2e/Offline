# Configuration file for runCosmic.py
#  - Generate N cosmic events (no other processes)
#  - Run these through G4.
#  - No event display.
#  - Write event data to an output file
#  - Save state of random numbers to the event-data output file
#
# $Id: runCosmic.py,v 1.2 2010/08/04 09:49:09 timothym Exp $
# $Author: timothym $
# $Date: 2010/08/04 09:49:09 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.  
process = mu2e.Process("runCosmic")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(1000)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("cosmics/cosmic.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberService",
                          globalSeed=mu2e.untracked.int32(9877)
))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_cosmic.txt")
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
process.generate = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/cosmicconfig_01.txt")
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    )

# Save state of random numbers to the event.
process.randomsaver = mu2e.EDAnalyzer("RandomNumberSaver")

# Define the output file.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:cosmics/data_cosmic.root'),
    outputCommands = cms.untracked.vstring(
     'keep *_*_*_*',
     'drop mu2eSimParticles_*_*_*'   # Uncomment this line to reduce file size.
    ),

)

# Look at the hits from G4.
process.checkhits = mu2e.EDFilter(
    "CosmicTuple",
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.001),
    traverseZ     = mu2e.double(1500.),
    minimump      = mu2e.double(103.),
    maximump      = mu2e.double(107.),
    pdgidnumber   = mu2e.double(11),
)

# End of the section that defines and configures modules.

# Adjust configuration of message logger.
# Enable debug printout from the module instance "hitinspect".
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.cerr.threshold = mu2e.untracked.string('DEBUG')
process.MessageLogger.debugModules.append("hitinspect")
process.MessageLogger.categories.append("ToyHitInfo")
process.MessageLogger.categories.append("GEOM")

# turn on timers
# process.Timing = mu2e.Service("Timing");

# turn on memory checker
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#   oncePerEventMode = cms.untracked.bool(False),
#   showMallocInfo = cms.untracked.bool(False),
#   ignoreTotal = cms.untracked.int32(100)
#)

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.generate*process.g4run*process.randomsaver*
                                process.checkhits*process.outfile );
