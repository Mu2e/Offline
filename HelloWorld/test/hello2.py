# Tell the framework to run the HelloWorld2 module.
#
# $Id: hello2.py,v 1.1 2010/09/01 18:56:46 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/09/01 18:56:46 $
#
# Original author Rob Kutschke
#
# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("HelloWorld2")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(2)
)

# Configure two modules.
process.source = mu2e.Source("EmptySource")
process.hello  = mu2e.EDAnalyzer("HelloWorld2",
       magicNumber = mu2e.untracked.int32(42)
)

# Tell the framework which modules to run.
process.output = mu2e.EndPath(process.hello);

