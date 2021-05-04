#ifndef Mu2eKinKal_KKFitSettings_hh
#define Mu2eKinKal_KKFitSettings_hh
//
// Struct for configuring the Mu2e KinKal fit 
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
#include "KinKal/Fit/Config.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
namespace mu2e {
  namespace Mu2eKinKal{

    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    // struct for defining the KinKal Config object
    struct KinKalConfig {
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level") };
      fhicl::Atom<int> minndof { Name("MinNDOF"), Comment("Minimum number of Degrees of Freedom to conitnue fitting") };
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration") };
      fhicl::Atom<float> dwt { Name("Deweight"), Comment("Deweighting factor when initializing the track end parameters") };
      fhicl::Atom<float> dparams { Name("DeltaParams"), Comment("Parameter difference threshold (units of chisquared)") };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for final fit (ns)") };
      fhicl::Atom<int> bfieldCorr { Name("BFieldCorrection"), Comment("BField correction algorithm") };
      fhicl::Atom<float> btol { Name("BCorrTolerance"), Comment("Tolerance on BField correction accuracy (mm)") };
      using MetaIterationSettings = fhicl::Sequence<fhicl::Tuple<float,float,float>>;
      MetaIterationSettings miConfig { Name("MetaIterationSettings"), Comment("MetaIteration sequence configuration parameters, format: \n"
      " 'Temperature (dimensionless)', 'Delta chisquared/DOF for convergence', 'Delta chisquared/DOF for divergence'") };
      using KKStrawHitUpdaterSettings = fhicl::Sequence<fhicl::Tuple<float,float,float,int,size_t>>;
      KKStrawHitUpdaterSettings shuConfig{ Name("KKStrawHitUpdaterSettings"), Comment("KKStrawHitUpdater settings, format: \n"
      " 'Minimum wire DOCA to assign L/R ambiguity and use drift'', 'Maximum wire DOCA to use hit', 'Maximum wire DOCA chito use hit', 'NullHit dimension','Meta-iteration'") };
    };
  // function to convert fhicl configuration to KinKal Config object
    KinKal::Config makeConfig(KinKalConfig const& fconfig);

  // struct for configuring KKFit object
    struct KKFitConfig {
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level") };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for final fit (ns)") };
      fhicl::Atom<float> tpocaPrec { Name("TPOCAPrecision"), Comment("TPOCA calculation precision (ns)") };
      fhicl::Atom<int> nullHitDimension { Name("NullHitDimension"), Comment("Null hit constrain dimension") }; 
      fhicl::Atom<float> nullHitVarianceScale { Name("NullHitVarianceScale"), Comment("Scale factor on geometric null hit variance") }; 
      fhicl::Atom<int> fitParticle {  Name("FitParticle"), Comment("Particle type to fit: e-, e+, mu-, ...")};
      fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream") };
      fhicl::Atom<bool> addMaterial { Name("AddMaterial"), Comment("Add material effects to the fit") }; 
      fhicl::Atom<bool> useCaloCluster { Name("UseCaloCluster"), Comment("Use CaloCluster in the fit") }; 
      fhicl::Atom<float> caloDt{ Name("CaloTrackerTimeOffset"), Comment("Time offset of calorimeter data WRT  (ns)") };
      fhicl::Atom<float> caloPosRes{ Name("CaloPositionResolution"), Comment("Transverse resolution of CaloCluster position (mm)") }; // this should come from the CaloCluster FIXME!
      fhicl::Atom<float> caloPropSpeed{ Name("CaloPropagationSpeed"), Comment("Axial speed of light in a crystal (mm/ns)") }; // see doc 25320.  this should come from the CaloCluster FIXME!
      fhicl::Atom<float> minCaloEnergy{ Name("MinCaloClusterEnergy"), Comment("Minimum CaloCluster energy to use as a KKCaloHit (MeV)") };
      fhicl::Atom<float> maxCaloDt{ Name("MaxCaloClusterDt"), Comment("Maximum CaloCluster time - track extrapolation time (ns)") };
      fhicl::Atom<float> maxCaloDoca { Name("MaxCaloClusterDOCA"), Comment("Max DOCA to add a CaloCluster (mm)") };
      fhicl::Sequence<std::string> addHitSelect { Name("AddHitSelect"), Comment("Flags required to be present to add a hit") };
      fhicl::Sequence<std::string> addHitReject { Name("AddHitReject"), Comment("Flags required not to be present to add a hit") };
      fhicl::Atom<float> maxStrawHitDOCA { Name("MaxStrawHitDOCA"), Comment("Max DOCA to add a hit (mm)") };
      fhicl::Atom<float> maxStrawHitDt { Name("MaxStrawHitDt"), Comment("Max Detla time to add a hit (ns)") };
      fhicl::Atom<float> maxStrawHitChi { Name("MaxStrawHitChi"), Comment("Max Chi to add a hit") };
      fhicl::Atom<int> strawBuffer { Name("StrawBuffer"), Comment("Buffer to add when searching for straws") };
      fhicl::Atom<float> maxStrawDOCA { Name("MaxStrawDOCA"), Comment("Max DOCA to add straw material (mm)") };
    };
  }
}
#endif
