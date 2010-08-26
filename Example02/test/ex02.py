# Configuration file for Example/02
#
# $Id: ex02.py,v 1.7 2010/08/26 15:50:25 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/08/26 15:50:25 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("Ex02")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(400)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("ex02histo.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Example02/test/ex02geom.txt")
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mu2e.Source("EmptySource")

# Make some hits and add them to the event.
process.ex02hitmaker = mu2e.EDProducer(
    "Ex02MakeHits",
    minPulseHeight=mu2e.double(6.5),
    seed=mu2e.untracked.vint32(7790)
    )

# Filter to select events with an odd event number.
process.oddfilter = mu2e.EDFilter(
    "Ex02SelectEvents",
    keepOddOrEven=mu2e.untracked.int32(1)
    )

# Filter to select events with an odd event number.
process.evenfilter = mu2e.EDFilter(
    "Ex02SelectEvents",
    keepOddOrEven=mu2e.untracked.int32(2)
    )

# Get the hits out of the event and look at them.  Make histograms.
process.hitinspect = mu2e.EDAnalyzer(
    "Ex02InspectHits",
    maxFullPrint=mu2e.untracked.int32(10)
    )

# Write the odd numbered events to one output file.
process.outputodd = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:ex02outputodd.root'),
    SelectEvents = mu2e.untracked.PSet(
        SelectEvents = mu2e.vstring('podd')
    )
)

# Write the even numbered event to a different output file.
process.outputeven = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:ex02outputeven.root'),
    SelectEvents = mu2e.untracked.PSet(
        SelectEvents = mu2e.vstring('peven')
    )
)

# End of the section that defines and configures modules.

# Define a path that selects odd numbered events
process.podd = mu2e.Path( process.ex02hitmaker
                         *process.oddfilter
                        )

# Define a path that selects even numbered events
process.peven = mu2e.Path( process.ex02hitmaker
                          *process.evenfilter
                         )

# Adjust configuration of message logger.
# Enable debug printout from the module instance "hitinspect".
# Print unlimited messages with category ToyHitInfo.
process.MessageLogger.cerr.threshold = mu2e.untracked.string('DEBUG')
process.MessageLogger.debugModules.append("hitinspect")
process.MessageLogger.categories.append("ToyHitInfo")
#process.MessageLogger.categories.append("GEOM")
process.MessageLogger.cerr.ToyHitInfo = mu2e.untracked.PSet(
    limit = mu2e.untracked.int32(-1)
)

# Tell the system to execute all paths.
process.output = mu2e.EndPath(   process.hitinspect
                               *(process.outputeven + process.outputodd)
                             )

