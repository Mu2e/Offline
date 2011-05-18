# Configuration file for Histforpabs
# Histforpabs the output of g4test_03.py; make histograms and printout.
#
# $Id: histforpabs.py,v 1.0 2010/09/16
# $Author: wb $

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("Histforpabs01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("histforpabs.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Initialize the random number sequences.
# This just changes the seed for the global CLHEP random engine.
process.add_(mu2e.Service("RandomNumberGeneratorService"))

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
    "Histforpabs",
    g4ModuleLabel = mu2e.string("g4run"),
    minimumEnergy = mu2e.double(0.001),
)

# End of the section that defines and configures modules.

# Tell the system to execute the module.
process.output = mu2e.EndPath(  process.checkhits );

