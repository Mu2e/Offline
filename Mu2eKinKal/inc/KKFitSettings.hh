#ifndef Mu2eKinKal_KKFitSettings_hh
#define Mu2eKinKal_KKFitSettings_hh
//
// Struct for configuring the Mu2e KinKal fit 
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
#include "DataProducts/inc/PDGCode.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
namespace mu2e {
  namespace Mu2eKinKal{
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
    struct FitSettings {
      using MetaIterationSettings = fhicl::Sequence<fhicl::Tuple<float,float,float>>;
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration"), 10 };
      fhicl::Atom<float> tpocaPrec { Name("TPOCAPrecision"), Comment("TPOCA calculation precision (ns)"), 1e-4 };
      fhicl::Atom<float> dwt { Name("Deweight"), Comment("Deweighting factor when initializing the track end parameters"), 1.0e6 };
      fhicl::Atom<float> dparams { Name("DeltaParams"), Comment("Parameter difference threshold (units of chisquared)"), 1.0e4 };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for final fit (ns)"), 0.1 };
      fhicl::Atom<float> btol { Name("BCorrTolerance"), Comment("Tolerance on BField correction accuracy (mm)"), 0.01 };
      fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };
      fhicl::Atom<int> bfieldCorr { Name("BFieldCorrection"), Comment("BField correction algorithm") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Atom<int> minndof { Name("MinNDOF"), Comment("Minimum number of Degrees of Freedom to conitnue fitting"), 5  };
      fhicl::Atom<int> nullHitDimension { Name("NullHitDimension"), Comment("Null hit constrain dimension"), 2 }; 
      fhicl::Atom<float> nullHitVarianceScale { Name("NullHitVarianceScale"), Comment("Scale factor on geometric null hit variance"), 1.0 }; 
      fhicl::Atom<int> fitParticle {  Name("FitParticle"), Comment("Particle type to fit: e-, e+, mu-, ..."), PDGCode::e_minus};
      fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream"), TrkFitDirection::downstream };
      fhicl::Atom<bool> addMaterial { Name("AddMaterial"), Comment("Add material effects to the fit"), true }; 
      fhicl::Atom<bool> useCaloCluster { Name("UseCaloCluster"), Comment("Use CaloCluster in the fit"), true }; 
      fhicl::Atom<float> caloDt{ Name("CaloTrackerTimeOffset"), Comment("Time offset of calorimeter data WRT tracker (ns)"), -0.1 };
      fhicl::Atom<float> caloPosRes{ Name("CaloPositionResolution"), Comment("Transverse resolution of CaloCluster position (mm)"), 15.0 };
      fhicl::Atom<float> caloPropSpeed{ Name("CaloPropagationSpeed"), Comment("Axial speed of light in a crystal (mm/ns)"), 200.0 }; // see doc 25320
      fhicl::Sequence<std::string> addHitSelect { Name("AddHitSelect"), Comment("Flags required to be present to add a hit"), std::vector<std::string>() };
      fhicl::Sequence<std::string> addHitReject { Name("AddHitReject"), Comment("Flags required not to be present to add a hit"), std::vector<std::string>() };
      fhicl::Atom<float> maxAddDOCA { Name("MaxAddDOCA"), Comment("Max DOCA to add a hit or straw (mm)"), 2.75 };
      fhicl::Atom<float> maxAddDt { Name("MaxAddDt"), Comment("Max Detla time to add a hit (ns)"), 3.0 };
      fhicl::Atom<float> maxAddChi { Name("MaxAddChi"), Comment("Max Chi to add a hit"), 4.0 };
      fhicl::Atom<float> maxAddDeltaU { Name("MaxAddDeltaU"), Comment("Max Delta-U to add a straw (mm)"), 10.0 };
      MetaIterationSettings mconfig { Name("MetaIterationSettings"), Comment("MetaIteration sequence configuration parameters, format: \n"
      " 'Temperature (dimensionless)', Delta chisquared/DOF for convergence', 'Delta chisquared/DOF for divergence'") };
    };
    // function to turn a StrawHit into a Line object
    class ComboHit;
    class Straw;
    class StrawResponse;
    KinKal::Line hitLine(ComboHit const& ch, Straw const& straw,StrawResponse const& strawresponse);
  }
}
#endif
