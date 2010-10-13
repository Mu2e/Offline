# Configuration file for G4Test04
# Similar to g4test_03.py but do fewer events and do not make the output file.
#
# $Id: g4test_04.py,v 1.8 2010/10/13 23:23:32 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/10/13 23:23:32 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.  
process = mu2e.Process("G4Test04")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(2)
)


# Define the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Define the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("g4test_04.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define the random number generator service.
process.RandomNumberGeneratorService = mu2e.Service("RandomNumberGeneratorService")

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
)

# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# Uncomment to enable per module timing
#process.Timing = mu2e.Service("Timing",
#    useJobReport = mu2e.untracked.bool(True)
#)

# Uncomment to enable memory use profiling
#process.SimpleMemoryCheck = mu2e.Service("SimpleMemoryCheck",
#    oncePerEventMode = mu2e.untracked.bool(False),
#    showMallocInfo = mu2e.untracked.bool(False),
#    ignoreTotal = mu2e.untracked.int32(5)
#)

# Uncomment to enable trace printout to show what the framework calls when.
# process.Tracer = mu2e.Service("Tracer")

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mu2e.Source("EmptySource")

#  Make some generated tracks and add them to the event.
process.generate = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/genconfig_02.txt"),
    seed=mu2e.untracked.vint32(7789)
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    seed=mu2e.untracked.vint32(9877)
    )

# Form StrawHits (SH).
process.makeSH = mu2e.EDProducer(
    "MakeStrawHit",
    g4ModuleLabel = mu2e.string("g4run"),
    seed=mu2e.untracked.vint32(7790),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Look at the hits from G4.
process.checkhits = mu2e.EDAnalyzer(
    "ReadBack",
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.001),
    maxFullPrint  = mu2e.untracked.int32(201)
)

# Save state of random numbers to the event.
process.randomsaver = mu2e.EDAnalyzer("RandomNumberSaver")

# End of the section that defines and configures modules.

# Adjust configuration of message logger.
# Enable debug printout from the module instance "hitinspect".
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.cerr.threshold = mu2e.untracked.string('DEBUG')
process.MessageLogger.debugModules.append("hitinspect")
process.MessageLogger.categories.append("ToyHitInfo")
process.MessageLogger.categories.append("GEOM")

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.generate*process.g4run*process.makeSH*process.randomsaver*
                                process.checkhits);
