#
# Run the StoppingTarget00 module.
#
# $Id: stoppingTarget00.py,v 1.2 2011/05/18 02:27:14 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:14 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("StoppingTarget00")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(20000)
)

# Define the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Define the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("stoppingTarget00.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/beamline_geom03a_readback.txt")
)

# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# Read events from a file
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("/prj/mu2e/users/kutschke/beamlineData_01.root")
)

#
process.stopping = mu2e.EDAnalyzer(
    "StoppingTarget00",
    g4ModuleLabel        = mu2e.string("g4run")
)

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.stopping )
