#
# Event display for a sample of conversion electrons.
#
# $Id: conversions.py,v 1.1 2011/03/10 20:41:17 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/03/10 20:41:17 $
#

import FWCore.ParameterSet.python.Config as mu2e

process = mu2e.Process("EvtDisplayConversions")

process.load("MessageLogger_cfi")

process.ConditionsService = mu2e.Service("ConditionsService",
        conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt"))

process.GeometryService = mu2e.Service("GeometryService",
        inputfile=mu2e.untracked.string("Mu2eG4/test/geom_cosmic.txt"))

process.source = mu2e.Source("PoolSource",
  fileNames = mu2e.untracked.vstring("/grid/fermiapp/mu2e/DataFiles/ExampleDataFiles/Tutorials/conversionOnly.root"))

process.eventdisplay = mu2e.EDAnalyzer("EventDisplay")

process.output = mu2e.EndPath(process.eventdisplay);

