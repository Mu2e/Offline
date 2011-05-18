import FWCore.ParameterSet.python.Config as mu2e

process = mu2e.Process("EventDisplay")

process.load("MessageLogger_cfi")

#Mu2eUtilities got changed and requires now loading ConditionsService,
#since the Mu2eUtilities shared library doesn't get linked with the ConditionsService library
process.ConditionsService = mu2e.Service("ConditionsService",
                            conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt"))

process.GeometryService = mu2e.Service("GeometryService",
                          inputfile=mu2e.untracked.string("Mu2eG4/test/geom_cosmic.txt"))

process.source = mu2e.Source("PoolSource",
                 fileNames = mu2e.untracked.vstring("/mu2e/data/outstage/wasiko/56545/56545_176/data_cosmic.root"))

process.eventdisplay = mu2e.EDAnalyzer("EventDisplay")

process.output = mu2e.EndPath(process.eventdisplay);

