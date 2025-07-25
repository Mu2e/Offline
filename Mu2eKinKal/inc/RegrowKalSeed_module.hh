//
// Recreate a KKTrack from a KalSeed. This can include additional extrapolation, new calibrations/alignment, etc.
// By default, no Pat. Rec. or straw hit state re-assignment is done.
//
// Original author: D. Brown (LBNL) 4/18/2025
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Tuple.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
// conditions
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/KinKalGeom/inc/SurfaceMap.hh"
// utiliites
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKConstantBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateToZ.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateIPA.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateST.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
// C++
#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <memory>

namespace mu2e {
  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;

  struct RegrowKalSeedConfig {
    fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedPtrCollection"), Comment("KalSeedPtr collection to processed ") };
    fhicl::Atom<art::InputTag> comboHitCollection {Name("ComboHitCollection"), Comment("Reduced ComboHit collection ") };
    fhicl::Atom<art::InputTag> indexMap {Name("StrawDigiIndexMap"), Comment("Map between original and reduced ComboHits") };
  };

  // Extrapolation configuration
  struct KKExtrapConfig {
    fhicl::Atom<float> Tol { Name("Tolerance"), Comment("Tolerance on fractional momemtum precision when extrapolating fits") };
    fhicl::Atom<float> MaxDt { Name("MaxDt"), Comment("Maximum time to extrapolate a fit") };
    fhicl::Sequence<int> SurfaceIDs { Name("Surfaces"), Comment("Surface IDs to extrapolate to") };
    fhicl::Sequence<int> ExtrapolationDirection { Name("ExtrapolationDirections"), Comment("Time Directions to extrapolate in") };
  };
  //
  template <class KTRAJ> class RegrowKalSeed : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<RegrowKalSeedConfig>;
      explicit RegrowKalSeed(const Parameters& settings);
  };
}
