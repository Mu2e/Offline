# First example of accessing Mu2e data.
#
# Read back the output of g4test_03.py; make a few histograms.
#
# $Id: readback0.py,v 1.2 2011/05/18 02:27:18 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:18 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("ReadBack00")

# Maximum number of events to do: -1 = read to end of file.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(-1)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("readback0.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define the geometry: file name should be the same as for generation.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/geom_01.txt")
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
    "ReadBack0",
    minimumEnergy     = mu2e.double(0.001),
)

# End of the section that defines and configures modules.

# Tell the system to execute the module.
process.output = mu2e.EndPath(  process.checkhits );

