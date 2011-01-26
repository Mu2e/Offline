import FWCore.ParameterSet.python.Config as mu2e

process = mu2e.Process("EventDisplay")

process.maxEvents = mu2e.untracked.PSet(input = mu2e.untracked.int32(100))

process.load("Config/MessageLogger_cfi")

process.GeometryService = mu2e.Service("GeometryService",
                          inputfile=mu2e.untracked.string("EventDisplay/test/geom_cosmic.txt"))

process.source = mu2e.Source("PoolSource",
                 fileNames = mu2e.untracked.vstring("EventDisplay/test/data_cosmic.root"))

process.root1 = mu2e.EDAnalyzer("EventDisplay", g4ModuleLabel = mu2e.string("g4run"))

process.output = mu2e.EndPath(process.root1);

