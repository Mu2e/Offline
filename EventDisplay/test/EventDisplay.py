import FWCore.ParameterSet.python.Config as mu2e

process = mu2e.Process("EventDisplay")

process.load("Config/MessageLogger_cfi")

process.GeometryService = mu2e.Service("GeometryService",
                          inputfile=mu2e.untracked.string("Mu2eG4/test/geom_cosmic.txt"))

process.source = mu2e.Source("PoolSource",
                 fileNames = mu2e.untracked.vstring("/mu2e/data/outstage/wasiko/56545/56545_8/data_cosmic.root"))

process.eventdisplay = mu2e.EDAnalyzer("EventDisplay")

process.output = mu2e.EndPath(process.eventdisplay);

