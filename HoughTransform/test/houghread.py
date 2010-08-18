# Configuration file for G4Test03
#
# $Id: houghread.py,v 1.2 2010/08/18 05:12:34 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/08/18 05:12:34 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("HoughRead0501")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(100)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("houghread_histos_501.root"),
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

# read from a file
process.source = mu2e.Source("PoolSource",
                 fileNames = mu2e.untracked.vstring("hough_0501.root") )

#  Make some generated tracks and add them to the event.
#1 is nothing but conversions
#2 has cosmics and radiative pi's

process.generate = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/genconfig_03.txt")
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("generate"),
    seed=mu2e.untracked.vint32(9877),
)


process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('file:houghread_0501.root'),
)

# apply noise and inefficiency
process.simpleen = mu2e.EDProducer(
    "SimpleEffyNoise",    
    diagLevel = mu2e.untracked.int32(1),    
    noiseRate = mu2e.double(.02),    
    hitIneffy = mu2e.double(.0) )

# Look at the hits from G4.
process.checkhits = mu2e.EDProducer(
# this line tells me where to get the plugin file from    
    "HoughTest",
    NPeakSearch = mu2e.uint32(50),
    hitCreatorName = mu2e.string("simpleen"),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Look at the hits from G4.
process.plothits = mu2e.EDAnalyzer(
# this line tells me where to get the plugin file from    
    "HoughTuner",
    hitCreatorName = mu2e.string("simpleen"),
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
process.output = mu2e.EndPath(  process.plothits );

