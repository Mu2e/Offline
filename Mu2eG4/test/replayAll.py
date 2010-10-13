# Configuration file to replay the output of g4test_03 using
# the random number state saved in the event.  The output
# event contains both sets of data products, those created in
# the first run and those created in the second run.
#
# $Id: replayAll.py,v 1.6 2010/10/13 23:38:47 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/10/13 23:38:47 $
#
# Original author Rob Kutschke
#
# Notes:
# 1) There is a bug in the interaction between root and the persistency
#    mechanism in the framework.  The problem is that fast cloning does
#    not work.  The solution is to disable fast cloning, which is done
#    in the configuration of the output modules.
#

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("replayAll")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("replayAll.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# State of random number engines will be restored from the input event.
# Define the random number generator service.
process.RandomNumberGeneratorService = mu2e.Service("RandomNumberGeneratorService",
                                       restoreStateLabel=mu2e.untracked.string("randomsaver")
)

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("drop_03.root"),
)

#  Make some generated tracks and add them to the event.
process.generate = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/genconfig_02.txt")
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    )

# Define the output file. See note 1.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:replayAllData.root'),
    fastCloning = mu2e.untracked.bool(False)
)

# Look at the hits from G4.
process.checkhits2 = mu2e.EDAnalyzer(
    "ReadBack",
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.001),
    maxFullPrint  = mu2e.untracked.int32(201)
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
process.output = mu2e.EndPath(  process.generate*process.g4run*process.checkhits2*process.outfile );
