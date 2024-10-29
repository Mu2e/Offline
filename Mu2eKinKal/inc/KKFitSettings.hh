#ifndef Mu2eKinKal_KKFitSettings_hh
#define Mu2eKinKal_KKFitSettings_hh
//
// Structs for configuring the Mu2e KinKal fit
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Tuple.h"
#include "KinKal/Fit/Config.hh"
#include "canvas/Utilities/InputTag.h"
#include "KinKal/Fit/MetaIterConfig.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
// updaters
#include "Offline/Mu2eKinKal/inc/CADSHU.hh"
#include "Offline/Mu2eKinKal/inc/DriftANNSHU.hh"
#include "Offline/Mu2eKinKal/inc/BkgANNSHU.hh"
#include "Offline/Mu2eKinKal/inc/Chi2SHU.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
namespace mu2e {
  namespace Mu2eKinKal{

    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    // struct for defining the KinKal Config object and updaters
    struct KinKalConfig {
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level") };
      fhicl::Atom<int> minndof { Name("MinNDOF"), Comment("Minimum number of Degrees of Freedom to conitnue fitting") };
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration") };
      fhicl::Atom<float> dwt { Name("Deweight"), Comment("Deweighting factor when initializing the track end parameters") };
      fhicl::Atom<float> convdchisq { Name("ConvergenceDeltaChisq"), Comment("Maximum Chisq/DOF change between iterations to define convergence") };
      fhicl::Atom<float> divdchisq { Name("DivergenceDeltaChisq"), Comment("Minimum Chisq/DOF change between iterations to define divergence") };
      fhicl::Atom<float> dparams { Name("DivergenceDeltaParams"), Comment("Parameter difference chisquared to define divergence threshold") };
      fhicl::Atom<float> dgap { Name("DivergenceGap"), Comment("Average trajectory gap to define divergence (mm)") };
      fhicl::Atom<bool> bfieldCorr { Name("BFieldCorrection"), Comment("Apply correction for BField inhomogeneity") };
      fhicl::Atom<bool> ends { Name("ProcessEnds"), Comment("Process purely passive sites at the time range ends") };
      fhicl::Atom<float> btol { Name("BCorrTolerance"), Comment("Tolerance on BField correction momentum fractional accuracy (dimensionless)") };
      // Updater settings
      using MetaIterationSettings = fhicl::Sequence<fhicl::Tuple<float,std::string>>;
      MetaIterationSettings miConfig { Name("MetaIterationSettings"), Comment("Temperature (dimensionless), StrawHitUpdater algorithm") };
      using CADSHUSettings = fhicl::OptionalSequence<fhicl::Tuple<float,float,float,float,std::string,std::string,std::string,int>>;
      CADSHUSettings cadshuConfig{ Name("CADSHUSettings"), Comment(CADSHU::configDescription()) };
      using DriftANNSHUSettings = fhicl::OptionalSequence<fhicl::Tuple<std::string,float,std::string,float,float,std::string,std::string,int>>;
      DriftANNSHUSettings annshuConfig{ Name("DriftANNSHUSettings"), Comment(DriftANNSHU::configDescription()) };
      using BkgANNSHUSettings = fhicl::OptionalSequence<fhicl::Tuple<std::string,float,std::string,int>>;
      BkgANNSHUSettings bkgshuConfig{ Name("BkgANNSHUSettings"), Comment(BkgANNSHU::configDescription()) };
      using Chi2SHUSettings = fhicl::OptionalSequence<fhicl::Tuple<unsigned,float,float,float,std::string,std::string,std::string,std::string,int>>;
      Chi2SHUSettings combishuConfig{ Name("Chi2SHUSettings"), Comment(Chi2SHU::configDescription()) };
      using StrawXingUpdaterSettings = fhicl::Sequence<fhicl::Tuple<float,float,float,bool,int>>;
      StrawXingUpdaterSettings sxuConfig{ Name("StrawXingUpdaterSettings"), Comment(StrawXingUpdater::configDescription()) };
    };
    // function to convert fhicl configuration to KinKal Config object
    KinKal::Config makeConfig(KinKalConfig const& fconfig);

    // struct for final fit configuration. This just changes convergence parameters
    struct KKFinalConfig {
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration") };
      fhicl::Atom<float> convdchisq { Name("ConvergenceDeltaChisq"), Comment("Maximum Chisq/DOF change between iterations to define convergence") };
    };

    // struct for configuring KKFit object
    struct KKFitConfig {
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level") };
      fhicl::Atom<float> tpocaPrec { Name("TPOCAPrecision"), Comment("TPOCA calculation precision (ns)") };
      fhicl::Atom<bool> matCorr { Name("MaterialCorrection"), Comment("Correct the fit fo material effects") };
      fhicl::Atom<bool> addHits { Name("AddHits"), Comment("Add hits to the fit") };
      fhicl::Atom<bool> addMaterial { Name("AddMaterial"), Comment("Add materials to the fit") };
      fhicl::Atom<bool> useCaloCluster { Name("UseCaloCluster"), Comment("Use CaloCluster in the fit") };
      fhicl::Atom<unsigned> minNStrawHits { Name("MinNStrawHits"), Comment("Minimum number of straw hits to attempt a fit") };
      fhicl::Atom<size_t> strawHitClusterDeltaStraw { Name("StrawHitClusterDeltaStraw"), Comment("Maximum straw index difference between StrawHits in StrawHitClusters") };
      fhicl::Atom<float> strawHitClusterDeltaT { Name("StrawHitClusterDeltaT"), Comment("Maximum time difference between StrawHits in StrawHitClusters") };
      fhicl::Atom<std::string> strawHitClusterLevel { Name("StrawHitClusterLevel"), Comment("Level for selecting StrawHitClusters (see StrawIdMask for details") };
      fhicl::Atom<float> caloDt{ Name("CaloTrackerTimeOffset"), Comment("Time offset of calorimeter data WRT tracker (ns)") }; // this should come from the database FIXME
      fhicl::Atom<float> caloPosRes{ Name("CaloPositionResolution"), Comment("Transverse resolution of CaloCluster position (mm)") }; // this should come from the CaloCluster FIXME!
      fhicl::Atom<float> caloTimeRes{ Name("CaloTimeResolution"), Comment("Effective resolution of CaloCluster time (ns)") }; // this should come from the CaloCluster FIXME!
      fhicl::Atom<float> caloPropSpeed{ Name("CaloPropagationSpeed"), Comment("Axial speed of light in a crystal (mm/ns)") }; // see doc 25320.  this should come from the CaloCluster FIXME!
      fhicl::Atom<float> minCaloEnergy{ Name("MinCaloClusterEnergy"), Comment("Minimum CaloCluster energy to use as a KKCaloHit (MeV)") };
      fhicl::Atom<float> maxCaloDt{ Name("MaxCaloClusterDt"), Comment("Maximum CaloCluster time - track extrapolation time (ns)") };
      fhicl::Atom<float> maxCaloDoca { Name("MaxCaloClusterDOCA"), Comment("Max DOCA to add a CaloCluster (mm)") };
      fhicl::Sequence<std::string> addHitSelect { Name("AddHitSelect"), Comment("Flags required to be present to add a hit") };
      fhicl::Sequence<std::string> addHitReject { Name("AddHitReject"), Comment("Flags required not to be present to add a hit") };
      fhicl::Atom<float> maxStrawHitDOCA { Name("MaxStrawHitDOCA"), Comment("Max DOCA to add a hit (mm)") };
      fhicl::Atom<float> maxStrawHitDt { Name("MaxStrawHitDt"), Comment("Max Detla time to add a hit (ns)") };
      fhicl::Atom<int> maxDStraw { Name("MaxDStraw"), Comment("Maximum (integer) straw separation when adding straw hits") };
      fhicl::Atom<float> maxStrawDOCA { Name("MaxStrawDOCA"), Comment("Max DOCA to add straw material (mm)") };
      fhicl::Atom<float> maxStrawDOCAConsistency { Name("MaxStrawDOCAConsistency"), Comment("Max DOCA chi-consistency to add straw material") };
      // extension and sampling
      fhicl::Atom<std::string> saveTraj { Name("SaveTrajectory"), Comment("How to save the trajectory in the KalSeed: None, Full, or T0 (just the t0 segment)") };
    };
    // struct for configuring a KinKal fit module
    struct KKModuleConfig {
      fhicl::Atom<int> fitParticle {  Name("FitParticle"), Comment("Particle type to fit: e-, e+, mu-, ...")};
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     caloClusterCollection     {Name("CaloClusterCollection"),     Comment("CaloCluster collection ") };
      fhicl::Sequence<std::string> seedFlags { Name("SeedFlags"), Comment("Flags required to be present to convert a seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
    };
  }
}
#endif
