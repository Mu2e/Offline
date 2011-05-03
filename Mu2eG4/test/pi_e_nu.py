# Configuration file for studying pi->e nu at rest as a calibration line.
#
# $Id: pi_e_nu.py,v 1.1 2011/05/03 05:06:38 kutschke Exp $
# $Author: kutschke $
# $Date: 2011/05/03 05:06:38 $
#
# Original author Rob Kutschke
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mu2e

# Give this process a name.
process = mu2e.Process("PiENu01")

# Maximum number of events to do.
process.maxEvents = mu2e.untracked.PSet(
    input = mu2e.untracked.int32(200)
)

# Define the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
process.load("MessageLogger_cfi")

# Define the service that manages root files for histograms.
process.TFileService = mu2e.Service("TFileService",
                       fileName = mu2e.string("pi_e_nu_01.root"),
                       closeFileFast = mu2e.untracked.bool(False)
)

# Define the random number generator service.
process.RandomNumberGeneratorService = mu2e.Service("RandomNumberGeneratorService")

# Define the geometry.
process.GeometryService = mu2e.Service("GeometryService",
       inputfile=mu2e.untracked.string("Mu2eG4/test/pi_e_nu_geom.txt")
)


# Access the conditions data.
process.ConditionsService = mu2e.Service("ConditionsService",
       conditionsfile=mu2e.untracked.string("Mu2eG4/test/conditions_01.txt")
)

# A helper for our interface with G4.
process.G4Helper = mu2e.Service("G4Helper")

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mu2e.Source("EmptySource")

#  Make some generated tracks and add them to the event.
process.readfile = mu2e.EDProducer(
    "EventGenerator",
    inputfile = mu2e.untracked.string("Mu2eG4/test/pi_e_nu_genconfig.txt"),
    seed=mu2e.untracked.vint32(7789)
)

process.decay = mu2e.EDProducer(
    "EplusFromStoppedPion",
    inputModuleLabel=mu2e.string("readfile"),
    czmin=mu2e.double(-1.),
    czmax=mu2e.double(+1.),
    seed=mu2e.untracked.vint32(1234),
#    doHistograms=mu2e.untracked.bool(1),
)

# Run G4 and add its hits to the event.
process.g4run = mu2e.EDProducer(
    "G4",
    generatorModuleLabel = mu2e.string("decay"),
    seed=mu2e.untracked.vint32(9877)
    )

# Form StrawHits (SH).
process.makeStrawHits = mu2e.EDProducer(
    "MakeStrawHit",
    g4ModuleLabel = mu2e.string("g4run"),
    seed=mu2e.untracked.vint32(7790),
    diagLevel    = mu2e.untracked.int32(0),
    maxFullPrint = mu2e.untracked.int32(5)
)

# Form CaloHits (APD hits)
process.CaloReadoutHitsMaker =  mu2e.EDProducer(
    "MakeCaloReadoutHits",
    diagLevel      = mu2e.untracked.int32(0),
    maxFullPrint   = mu2e.untracked.int32(201),
    g4ModuleLabel = mu2e.string("g4run")
)

# Form CaloCrystalHits (reconstruct crystals from APDs)
process.CaloCrystalHitsMaker =  mu2e.EDProducer(
    "MakeCaloCrystalHits",
    diagLevel      = mu2e.untracked.int32(0),
    maxFullPrint   = mu2e.untracked.int32(201),
    g4ModuleLabel  = mu2e.string("g4run"),
    minimumEnergy  = mu2e.untracked.double(0.0),
    maximumEnergy  = mu2e.untracked.double(1000.0),
    minimumTimeGap = mu2e.untracked.double(100.0)
)

# Look at the hits from G4.
process.checkhits = mu2e.EDAnalyzer(
    "ReadBack",
    diagLevel            = mu2e.untracked.int32(0),
    g4ModuleLabel        = mu2e.string("g4run"),
    generatorModuleLabel = mu2e.string("decay"),
    minimumEnergy        = mu2e.double(0.001),
    maxFullPrint         = mu2e.untracked.int32(201)
)

process.readStrawHits = mu2e.EDAnalyzer("ReadStrawHit",
    makerModuleLabel = mu2e.string("makeStrawHits"),
    diagLevel    = mu2e.untracked.int32(3),
    maxFullPrint = mu2e.untracked.int32(100)
)


# Save state of random numbers to the event.
process.randomsaver = mu2e.EDAnalyzer("RandomNumberSaver")

# Define the output file.
process.outfile = mu2e.OutputModule(
    "PoolOutputModule",
    fileName = mu2e.untracked.string('file:data_pi_e_nu_01.root'),
    outputCommands = mu2e.untracked.vstring(
     'keep *_*_*_*',
     'drop mu2ePointTrajectoryMapVector_*_*_*',
#     'drop mu2eSimParticles_*_*_*'   # Uncomment this line to reduce file size.
    ),

)

# End of the section that defines and configures modules.

# Tell the system to execute all paths.
process.output = mu2e.EndPath(  process.readfile*
                                process.decay*
                                process.g4run*
                                process.makeStrawHits*
                                process.CaloReadoutHitsMaker*
                                process.CaloCrystalHitsMaker*
                                process.randomsaver*
                                process.checkhits*
                                process.readStrawHits*
                                process.outfile );
