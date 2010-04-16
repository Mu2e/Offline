# Tell the framework to run the hello world module.
#
# $Id: hello.py,v 1.1 2010/04/16 15:13:00 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/04/16 15:13:00 $
#
# Original author Rob Kutschke
#
# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("HelloWorld")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(2)
)

# Configure the modules.
process.source = mu2e.Source("EmptySource")
process.hello  = mu2e.EDAnalyzer("HelloWorld",
       magicNumber = mu2e.untracked.int32(42)
)

# The the framework which modules to run.
process.output = mu2e.EndPath(process.hello);

