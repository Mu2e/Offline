//
// Instantiation of RegrowKalSeed for LoopHelix fits
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
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KinKal::DVEC;
  using KinKal::TimeDir;
  using MatEnv::DetMaterial;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eKinKal::KKFinalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;

  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;
  struct RegrowLoopHelixConfig {
    fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedPtrCollection"), Comment("KalSeedPtr collection to processed ") };
    fhicl::Atom<art::InputTag> comboHitCollection {Name("ComboHitCollection"), Comment("Reduced ComboHit collection ") };
    fhicl::Atom<art::InputTag> indexMap {Name("StrawDigiIndexMap"), Comment("Map between original and reduced ComboHits") };

    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("RefitSettings") };
//     fhicl::OptionalTable<KKExtrapConfig> Extrapolation { Name("Extrapolation") }; needs to be pulled out of LoopHelixFit if needed
  };

  class RegrowLoopHelix : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<RegrowLoopHelixConfig>;
      using KTRAJ = KinKal::LoopHelix;
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KKTRK = KKTrack<KTRAJ>;
      using KKTRKCOL = OwningPointerCollection<KKTRK>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KKSTRAWXING = KKStrawXing<KTRAJ>;
      using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
      using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
      using KKIPAXING = KKShellXing<KTRAJ,KinKal::Cylinder>;
      using KKIPAXINGPTR = std::shared_ptr<KKIPAXING>;
      using KKIPAXINGCOL = std::vector<KKIPAXINGPTR>;
      using KKSTXING = KKShellXing<KTRAJ,KinKal::Annulus>;
      using KKSTXINGPTR = std::shared_ptr<KKSTXING>;
      using KKSTXINGCOL = std::vector<KKSTXINGPTR>;
      using KKCALOHIT = KKCaloHit<KTRAJ>;
      using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
      using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
      using KKFIT = KKFit<KTRAJ>;

      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;

      using KKMaterialConfig = KKMaterial::Config;

      explicit RegrowLoopHelix(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
      void endJob() override;
    private:
  };
  RegrowLoopHelix::RegrowLoopHelix(const Parameters& settings) : art::EDProducer(settings) {
    produces<KKTRKCOL>();
    produces<KalSeedCollection>();
  }

  void RegrowLoopHelix::beginRun(art::Run& run)
  {}

  void RegrowLoopHelix::produce(art::Event& event)
  {}

  void RegrowLoopHelix::endJob()
  {}

}
DEFINE_ART_MODULE(mu2e::RegrowLoopHelix)
