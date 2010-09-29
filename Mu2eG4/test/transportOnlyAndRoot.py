# Variant of g4test_03 but with transport only.
#
# $Id: transportOnlyAndRoot.py,v 1.1 2010/09/29 02:56:21 genser Exp $
# $Author: genser $
# $Date: 2010/09/29 02:56:21 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.  
process = mu2e.Process("transportAndRoot")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("transportOnlyAndRoot.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/transportOnlyGeom.txt")
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
    inputfile = mu2e.untracked.string("Mu2eG4/test/genconfig_tonly.txt"),
    seed=mu2e.untracked.vint32(7789)
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    rmvlevel = mu2e.untracked.int32(2),
#    visMacro = mu2e.untracked.string("Mu2eG4/test/vis45.mac"),
    seed=mu2e.untracked.vint32(9877)
)

# Save state of random numbers to the event.
process.randomsaver = mu2e.EDAnalyzer("RandomNumberSaver")

# Define the output file.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:data_06.root'),
    outputCommands = cms.untracked.vstring(
     'keep *_*_*_*',
#     'drop mu2eSimParticles_*_*_*'   # Uncomment this line to reduce file size.
    ),

)

# Look at the hits from G4.
process.checkhits = mu2e.EDAnalyzer(
    "ReadBack",
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.0),
    maxFullPrint  = mu2e.untracked.int32(201)
)


# Look at the geometry
process.geometryplots = mu2e.EDAnalyzer(
    "TTrackerGeomIntRootPlots",
)

# process.Tracer = mu2e.Service("Tracer")


# End of the section that defines and configures modules.

# Adjust configuration of message logger.
# Enable debug printout from the module instance "hitinspect".
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.cerr.threshold = mu2e.untracked.string('DEBUG')
process.MessageLogger.debugModules.append("hitinspect")
process.MessageLogger.categories.append("ToyHitInfo")
process.MessageLogger.categories.append("GEOM")

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.generate*process.g4run*process.randomsaver*process.geometryplots*
                                process.checkhits*process.outfile );
