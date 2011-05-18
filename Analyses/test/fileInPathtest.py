#
# Test of FileInPath.
#
# $Id: fileInPathtest.py,v 1.2 2011/05/18 02:27:14 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:14 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this job a name.
process = mu2e.Process("FileInPathTester")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(1)
)

process.source = mu2e.Source("EmptySource")

process.fiptest = mu2e.EDAnalyzer(
    "FileInPathTest",
    inputfile = mu2e.FileInPath("Mu2eG4/test/geom_01.txt"),
)

# End of the section that defines and configures modules.

# Tell the system to execute the modules.
process.output = mu2e.EndPath(  process.fiptest );
