// Configuration parameters used by Mu2eG4 and Mu2eG4MT modules
//
// Andrei Gaponenko, 2020


#ifndef Mu2eG4_inc_Mu2eG4Config_hh
#define Mu2eG4_inc_Mu2eG4Config_hh

#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTuple.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {
  namespace Mu2eG4Config {

    struct Debug {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> warnEveryNewRun {Name("warnEveryNewRun"), false };
      fhicl::Atom<bool> exportPDTStart {Name("exportPDTStart"), false };
      fhicl::Atom<bool> exportPDTEnd {Name("exportPDTEnd"), false };
      fhicl::Atom<std::string> storePhysicsTablesDir {Name("storePhysicsTablesDir"), ""};
      fhicl::Atom<int> diagLevel {Name("diagLevel"), 0};
      fhicl::Atom<int> worldVerbosityLevel {Name("worldVerbosityLevel"), 0};
      fhicl::Atom<int> steppingVerbosityLevel {Name("steppingVerbosityLevel"), 0};
      fhicl::Atom<int> trackingVerbosityLevel {Name("trackingVerbosityLevel"), 0};
      fhicl::Atom<bool> navigatorCheckMode {Name("navigatorCheckMode"), false};
      fhicl::Atom<int> navigatorVerbosityLevel {Name("navigatorVerbosityLevel"), 0};
      fhicl::Atom<int> PiENuPolicyVerbosity {Name("PiENuPolicyVerbosity"), 0};
      fhicl::Atom<int> mtDebugOutput {Name("mtDebugOutput"), 0};

      // fhicl::OptionalTuple<int,double,int> ionToGenerate { Name("ionToGenerate"), std::tuple<int,double>{1000591349, 0.163100 , 1} }; // particle code, excitation energy float level base index; experts only
      fhicl::OptionalTuple<int,double,int> ionToGenerate { Name("ionToGenerate") };

      fhicl::Atom<int> checkFieldMap {Name("checkFieldMap"), 0 };

      fhicl::Atom<bool> printElements {Name("printElements"), Comment("Print elements from constructMaterials()")};
      fhicl::Atom<bool> printMaterials {Name("printMaterials"), Comment("Print materials from constructMaterials()")};

      fhicl::Atom<bool> writeGDML {Name("writeGDML")};
      fhicl::Atom<std::string> GDMLFileName {Name("GDMLFileName")};

      fhicl::Atom<bool> stepLimitKillerVerbose {Name("stepLimitKillerVerbose")};
      fhicl::Sequence<int> eventList {Name("eventList"), std::vector<int>()};
      fhicl::Sequence<int> trackList {Name("trackList"), std::vector<int>()};
      fhicl::Sequence<int> trackingActionEventList {Name("trackingActionEventList"), std::vector<int>()};
      fhicl::Atom<bool> printTrackTiming {Name("printTrackTiming")};

    };

    struct Visualization {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> initMacro {Name("initMacro")};
      fhicl::Atom<std::string> GUIMacro {Name("GUIMacro")};
    };

    struct TimeVD_ {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Sequence<double> times {Name("times")};
      bool enabled() const { return !times().empty(); }
    };

    struct SDConfig_ {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Table<TimeVD_> TimeVD {Name("TimeVD")};

      fhicl::Atom<int> verbosityLevel {Name("verbosityLevel"), 0};

      fhicl::OptionalAtom<bool> enableAllSDs {Name("enableAllSDs"),
          Comment("Either set enableAllSDs=true or provide an explicit list in enableSD.")
          };

      fhicl::OptionalSequence<std::string> enableSD {Name("enableSD"),
          Comment("An explicit list of enabled detectors.  Must not be set if enableAllSDs=true.")
      };

      fhicl::Sequence<std::string> sensitiveVolumes {Name("sensitiveVolumes"), {}};
      fhicl::Sequence<std::string> preSimulatedHits {Name("preSimulatedHits"), {}};

      // FIXME: why is this necessary?
      fhicl::Sequence<std::string> inputs {Name("inputs"), {}};
      fhicl::Atom<double> cutMomentumMin {Name("cutMomentumMin"), 0.};
      fhicl::Atom<size_t> minTrackerStepPoints {Name("minTrackerStepPoints"), 15};
    };

    struct Physics {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      using OptionalDelegatedParameter = fhicl::OptionalDelegatedParameter;

      fhicl::Atom<std::string> stepper {Name("stepper")};
      fhicl::Atom<double> epsilonMin {Name("epsilonMin")};
      fhicl::Atom<double> epsilonMax {Name("epsilonMax")};
      fhicl::Atom<double> deltaOneStep {Name("deltaOneStep"), Comment("In mm")};
      fhicl::Atom<double> deltaIntersection {Name("deltaIntersection"), Comment("In mm")};
      fhicl::Atom<double> deltaChord {Name("deltaChord"), Comment("In mm")};
      fhicl::Atom<double> stepMinimum {Name("stepMinimum"), Comment("In mm")};
      fhicl::Atom<int> maxIntSteps {Name("maxIntSteps")};
      fhicl::Atom<double> bfieldMaxStep {Name("bfieldMaxStep"), Comment("In mm")};
      fhicl::Atom<double> strawGasMaxStep {Name("strawGasMaxStep"), Comment("In mm")};
      fhicl::Atom<bool> limitStepInAllVolumes {Name("limitStepInAllVolumes")};
      fhicl::Atom<bool> useEmOption4InTracker {Name("useEmOption4InTracker"), false};

      fhicl::Atom<double> protonProductionCut {Name("protonProductionCut")};

      fhicl::Atom<std::string> physicsListName {Name("physicsListName")};
      fhicl::Atom<bool> turnOffRadioactiveDecay {Name("turnOffRadioactiveDecay"),
          Comment("turnOn/Off radioactive decay. If both are false, we get the G4 default,\nwhich may be physics list dependent.\nIf both are true, it is a configuration error."), false};
      fhicl::Atom<bool> turnOnRadioactiveDecay {Name("turnOnRadioactiveDecay"), Comment("See turnOnRadioactiveDecay (above)"),false};
      fhicl::Atom<bool> turnOnThermalNeutronPhysics {Name("turnOnThermalNeutronPhysics"), false};
      fhicl::Atom<bool> decayMuonsWithSpin {Name("decayMuonsWithSpin"), false};
      fhicl::Atom<double> minRangeCut {Name("minRangeCut")};


      fhicl::Atom<bool> setMuHadLateralDisplacement {Name("setMuHadLateralDisplacement"), false};

      fhicl::Atom<bool> useDensityEffectInIonizationLossCalc {Name("useDensityEffectInIonizationLossCalc"), false};
      fhicl::Sequence<std::string> conductingMaterials {Name("conductingMaterials"),
          Comment("List of materials which are electrical conductors.")};

      fhicl::Sequence<int> noDecay {Name("noDecay"), Comment("List of PDG IDs that for which to turn decays off.")};

      fhicl::Atom<std::string> captureDModel {Name("captureDModel")};

      fhicl::Sequence<std::string> addProcesses {Name("addProcesses"), Comment("List process names.")};

      fhicl::Atom<std::string> PiENuPolicy {Name("PiENuPolicy"),
          Comment("pienu branching.  A string convertable to a double, "
                  "or one of \"PDG\", \"All\", \"None\" pre-defined settings.")
          };

      fhicl::OptionalAtom<bool> disableEnergyLossFluctuations {Name("disableEnergyLossFluctuations")};
#if G4VERSION>4110
      fhicl::OptionalAtom<unsigned int> setEnergyLossFluctuationModel {Name("energyLossFluctuationModel")};
#endif
      fhicl::OptionalAtom<double> mscModelTransitionEnergy {Name("mscModelTransitionEnergy")};
      fhicl::OptionalAtom<double> muonPreAssignedDecayProperTime {Name("muonPreAssignedDecayProperTime")};
      fhicl::OptionalAtom<double> muonMaxPreAssignedDecayProperTime {Name("muonMaxPreAssignedDecayProperTime")};
      fhicl::OptionalAtom<double> muonMinPreAssignedDecayProperTime {Name("muonMinPreAssignedDecayProperTime")};

#if G4VERSION>4112
      // works for Geant4 11.3.p02
      fhicl::OptionalAtom<bool> setBertiniAs11_2 {Name("setBertiniAs11_2")};
#endif

      OptionalDelegatedParameter BirksConsts {Name("BirksConsts")};
      OptionalDelegatedParameter minRangeRegionCuts {Name("minRangeRegionCuts")};

      fhicl::Atom<double> rangeToIgnore {Name("rangeToIgnore")};
    };

    struct Limits {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<unsigned> maxStepsPerTrack {Name("maxStepsPerTrack")};
      fhicl::Atom<unsigned> maxStepPointCollectionSize {Name("maxStepPointCollectionSize")};
      fhicl::Atom<unsigned> maxSimParticleCollectionSize {Name("maxSimParticleCollectionSize")};
    };

    struct TrajectoryControl_ {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<bool> produce {Name("produce")};
      fhicl::OptionalAtom<double> defaultMinPointDistance {Name("defaultMinPointDistance")};
      fhicl::OptionalAtom<unsigned> mcTrajectoryMinSteps {Name("mcTrajectoryMinSteps")};
      fhicl::OptionalAtom<double> mcTrajectoryMomentumCut {Name("mcTrajectoryMomentumCut")};
      fhicl::OptionalAtom<double> saveTrajectoryMomentumCut {Name("saveTrajectoryMomentumCut")};

      fhicl::OptionalDelegatedParameter perVolumeMinDistance {Name("perVolumeMinDistance"),
          Comment("A table that maps names to min distance between saved trajectory points.")
          };
    };

    struct EventLevelVolInfos {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> input { Name("input"), Comment("Tag of the event-level PhysicalVolumeInfoMultiCollection to read")};
      fhicl::Atom<std::string> outInstance { Name("outInstance"),
          Comment("Instance name of PhysicalVolumeInfoMultiCollection written into events."),
          "eventlevel"
          };
    };

    struct Inputs_ {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<std::string> primaryType { Name("primaryType"), Comment(primaryType_docstring()) };
      static std::string primaryType_docstring();
      bool multiStage() const { return primaryType() != "GenParticles"; }

      fhicl::Atom<art::InputTag> primaryTag {Name("primaryTag"), Comment("Tag of the input collection for G4 primaries")};

      fhicl::Atom<art::InputTag> inputMCTrajectories {Name("inputMCTrajectories"),
          Comment("MCTrajectoryCollection from the previous simulation stage, required for non-GenParticles primaries"),
          fhicl::use_if(this, &Inputs_::multiStage)
          };

      fhicl::OptionalAtom<unsigned> simStageOverride {Name("simStageOverride"),
          Comment("Normally simStage is determined at begin SubRun from the inputPhysVolumeMultiInfo product.\n"
                  "If simStageOverride is defined it will be used instead and inputPhysVolumeMultiInfo will only\n"
                  "be retrieved at end SubRun.  This option can only be used for non-GenParticles primaries."
                  ),
          fhicl::use_if(this, &Inputs_::multiStage)
          };

      fhicl::Atom<art::InputTag> inputPhysVolumeMultiInfo {Name("inputPhysVolumeMultiInfo"),
          Comment("phys volumes from the previous simulation stage, required for non-GenParticles primaries"),
          fhicl::use_if(this, &Inputs_::multiStage)
          };

      fhicl::OptionalTable<EventLevelVolInfos> updateEventLevelVolumeInfos {Name("updateEventLevelVolumeInfos"),
          Comment("An option to update event-level volume infos, for mixing jobs.\n"
                  "May only be used for non-GenParticles primaries."),
          fhicl::use_if(this, &Inputs_::multiStage)
          };

    };

    struct Top {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      using DelegatedParameter = fhicl::DelegatedParameter;

      fhicl::Table<Inputs_> inputs { Name("inputs") };

      fhicl::Table<Limits> ResourceLimits { Name("ResourceLimits") };
      fhicl::Table<TrajectoryControl_> TrajectoryControl { Name("TrajectoryControl") };
      fhicl::Table<SDConfig_> SDConfig { Name("SDConfig") };

      DelegatedParameter  Mu2eG4StackingOnlyCut { Name("Mu2eG4StackingOnlyCut") };
      DelegatedParameter  Mu2eG4SteppingOnlyCut { Name("Mu2eG4SteppingOnlyCut") };
      DelegatedParameter  Mu2eG4CommonCut { Name("Mu2eG4CommonCut") };

      fhicl::OptionalTable<SimParticleCollectionPrinter::Config> SimParticlePrinter { Name("SimParticlePrinter") };

      fhicl::Table<Physics> physics { Name("physics") };
      fhicl::Table<Debug> debug { Name("debug") };

      fhicl::Table<Visualization>  visualization { Name("visualization") };
      fhicl::Atom<std::string> g4Macro {Name("g4Macro"), ""};

      fhicl::Atom<std::string> generatorModuleLabel {Name("generatorModuleLabel"), ""};

      fhicl::Atom<std::string> salt {Name("salt"), ""};

      fhicl::Atom<bool> G4InteralFiltering {Name("G4InteralFiltering"), false};
    };
  }
}

#endif /*Mu2eG4_inc_Mu2eG4Config_hh*/
