#
# Print information about all data products in a file.
#
# $Id: dumpDataProducts.py,v 1.1 2011/02/23 07:03:22 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/02/23 07:03:22 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.  
process = mu2e.Process("DumpDataProducts")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(1)
)

# Load the standard message logger configuration.
process.load("MessageLogger_cfi")

# Read events from a file.
process.source = mu2e.Source("PoolSource",
   fileNames = mu2e.untracked.vstring("data_03.root")
)

# Print info about all data products in the file.
process.dump = mu2e.EDAnalyzer(
    "DataProductDump"
)

process.output = mu2e.EndPath(  process.dump );

