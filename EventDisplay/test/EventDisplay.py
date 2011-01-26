import FWCore.ParameterSet.python.Config as mu2e

process = mu2e.Process("EventDisplay")

process.maxEvents = mu2e.untracked.PSet(input = mu2e.untracked.int32(100))

process.load("Config/MessageLogger_cfi")

process.GeometryService = mu2e.Service("GeometryService",
                          inputfile=mu2e.untracked.string("Mu2eG4/test/geom_cosmic.txt"))

process.source = mu2e.Source("PoolSource",
                 fileNames = mu2e.untracked.vstring("/mu2e/data/outstage/wasiko/56545/56545_1/data_cosmic.root"))

process.root1 = mu2e.EDAnalyzer("EventDisplay", g4ModuleLabel = mu2e.string("g4run"))

process.output = mu2e.EndPath(process.root1);

