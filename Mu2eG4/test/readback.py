# Configuration file for Readback
# Readback the output of g4test_03.py; make histograms and printout.
#
# $Id: readback.py,v 1.11 2011/05/03 03:00:58 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/03 03:00:58 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("ReadBack01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("readback.root"),
                       closeFileFast = mu2e.untracked.bool(False)
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

# Read events from a file (made by example 3)
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Look at the hits from G4.
#  - minimum energy is in MeV
process.checkhits = mu2e.EDAnalyzer(
    "ReadBack",
    generatorModuleLabel=mu2e.string("generate"),
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.001),
    maxFullPrint = mu2e.untracked.int32(201)
)

# End of the section that defines and configures modules.

# Tell the system to execute the module.
process.output = mu2e.EndPath(  process.checkhits );

