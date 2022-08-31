#ifndef Mu2eKinKal_KKFitSettings_hh
#define Mu2eKinKal_KKFitSettings_hh
//
// Structs for configuring the Mu2e KinKal fit
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
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
      using MetaIterationSettings = fhicl::Sequence<fhicl::Tuple<float,int>>;
      MetaIterationSettings miConfig { Name("MetaIterationSettings"), Comment("MetaIteration sequence configuration parameters, format: \n"
          " Temperature (dimensionless), StrawHitUpdater algorithm") };
      using CAStrawHitUpdaterSettings = fhicl::OptionalSequence<fhicl::Tuple<float,float,float,float,std::string>>;
      CAStrawHitUpdaterSettings cashuConfig{ Name("CAStrawHitUpdaterSettings"), Comment("CAStrawHitUpdater settings, format: \n"
          " Minimum DOCA to use L/R ambiguity, Maximum wire DOCA to use hit, Minimum Dt to use hit, Maximum Dt to use a hit, States to freeze") };
      using ANNStrawHitUpdaterSettings = fhicl::OptionalSequence<fhicl::Tuple<std::string,float,float,std::string>>;
      ANNStrawHitUpdaterSettings annshuConfig{ Name("ANNStrawHitUpdaterSettings"), Comment("ANNStrawHitUpdater settings, format: \n"
          " Weight file, ann cut,  hit mode, null hit doca, states to freeze") };
      using BkgStrawHitUpdaterSettings = fhicl::OptionalSequence<fhicl::Tuple<std::string,float,std::string>>;
      BkgStrawHitUpdaterSettings bkgshuConfig{ Name("BkgStrawHitUpdaterSettings"), Comment("BkgStrawHitUpdater settings, format: \n"
          " Weight file, ann cut, states to freeze") };
      using CombinatoricStrawHitUpdaterSettings = fhicl::OptionalSequence<fhicl::Tuple<unsigned,float,float,float,float,std::string,std::string,int>>;
      CombinatoricStrawHitUpdaterSettings chuConfig{ Name("CombinatoricStrawHitUpdaterSettings"), Comment("CombinatoricStrawHitUpdater settings, format: \n"
          "Min Cluster Size, Inactive hit x^2 penalty, Null ambiguity x^2 penalty, Minimum significant x^2 difference, minimum drift DOCA, allowed states, states to freeze, diag level") };
      using StrawXingUpdaterSettings = fhicl::Sequence<fhicl::Tuple<float,float,float,bool>>;
      StrawXingUpdaterSettings sxuConfig{ Name("StrawXingUpdaterSettings"), Comment("StrawXingUpdater settings, format: \n"
          " Maximum DOCA to use straw material, Maximum DOCA to use unaveraged material,Maximum DOCA error to use unaveraged material, scale variance with annealing temp?") };
      //NB: when new updaters are introduced their config must be added as a new tuple sequences
    };
    // function to convert fhicl configuration to KinKal Config object
    KinKal::Config makeConfig(KinKalConfig const& fconfig);

    // struct for configuring Mu2e-specific fit settings
    struct Mu2eConfig {
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level") };
      fhicl::Atom<float> tpocaPrec { Name("TPOCAPrecision"), Comment("TPOCA calculation precision (ns)") };
      fhicl::Atom<int> fitParticle {  Name("FitParticle"), Comment("Particle type to fit: e-, e+, mu-, ...")};
      fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream") };
      fhicl::Atom<bool> addMaterial { Name("AddMaterial"), Comment("Add material effects to the fit") };
      fhicl::Atom<bool> useCaloCluster { Name("UseCaloCluster"), Comment("Use CaloCluster in the fit") };
      fhicl::Atom<size_t> strawHitClusterDeltaStraw { Name("StrawHitClusterDeltaStraw"), Comment("Maximum straw index difference between StrawHits in StrawHitClusters") };
      fhicl::Atom<float> strawHitClusterDeltaT { Name("StrawHitClusterDeltaT"), Comment("Maximum time difference between StrawHits in StrawHitClusters") };
      fhicl::Atom<std::string> strawHitClusterLevel { Name("StrawHitClusterLevel"), Comment("Level for selecting StrawHitClusters (see StrawIdMask for details") };
      fhicl::Atom<float> caloDt{ Name("CaloTrackerTimeOffset"), Comment("Time offset of calorimeter data WRT tracker (ns)") }; // this should come from the database FIXME
      fhicl::Atom<float> caloPosRes{ Name("CaloPositionResolution"), Comment("Transverse resolution of CaloCluster position (mm)") }; // this should come from the CaloCluster FIXME!
      fhicl::Atom<float> caloPropSpeed{ Name("CaloPropagationSpeed"), Comment("Axial speed of light in a crystal (mm/ns)") }; // see doc 25320.  this should come from the CaloCluster FIXME!
      fhicl::Atom<float> minCaloEnergy{ Name("MinCaloClusterEnergy"), Comment("Minimum CaloCluster energy to use as a KKCaloHit (MeV)") };
      fhicl::Atom<float> maxCaloDt{ Name("MaxCaloClusterDt"), Comment("Maximum CaloCluster time - track extrapolation time (ns)") };
      fhicl::Atom<float> maxCaloDoca { Name("MaxCaloClusterDOCA"), Comment("Max DOCA to add a CaloCluster (mm)") };
      fhicl::Sequence<std::string> addHitSelect { Name("AddHitSelect"), Comment("Flags required to be present to add a hit") };
      fhicl::Sequence<std::string> addHitReject { Name("AddHitReject"), Comment("Flags required not to be present to add a hit") };
      fhicl::Atom<float> maxStrawHitDOCA { Name("MaxStrawHitDOCA"), Comment("Max DOCA to add a hit (mm)") };
      fhicl::Atom<float> maxStrawHitDt { Name("MaxStrawHitDt"), Comment("Max Detla time to add a hit (ns)") };
      fhicl::Atom<int> strawBuffer { Name("StrawBuffer"), Comment("Buffer to add when searching for straws") };
      fhicl::Atom<float> maxStrawDOCA { Name("MaxStrawDOCA"), Comment("Max DOCA to add straw material (mm)") };
   };
  }
}
#endif
