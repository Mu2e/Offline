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
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/DataProducts/inc/IndexMap.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
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
#include "Offline/Mu2eKinKal/inc/KKExtrap.hh"
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
  using KinKal::Config;
  using MatEnv::DetMaterial;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eKinKal::KKFinalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;
  using KKMaterialConfig = KKMaterial::Config;
  using SDIS = std::set<StrawDigiIndex>;

  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;
  struct RegrowLoopHelixConfig {
    fhicl::Atom<int> debug{Name("debug"), Comment("Debug level"), 0};
    fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedCollection"), Comment("KalSeed collection to process") };
    fhicl::Atom<art::InputTag> comboHitCollection {Name("ComboHitCollection"), Comment("Reduced ComboHit collection") };
    fhicl::Atom<art::InputTag> caloClusterCollection {Name("CaloClusterCollection"),     Comment("CaloCluster collection ") };
    fhicl::Atom<art::InputTag> indexMap {Name("StrawDigiIndexMap"), Comment("Map between original and reduced ComboHits") };
    fhicl::OptionalAtom<art::InputTag> kalSeedMCAssns {Name("KalSeedMCAssns"), Comment("Association to KalSeedMC. If set, regrown KalSeeds will be associated with the same KalSeedMC as the original") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("RefitSettings") };
    fhicl::Atom<bool> extend {Name("Extend"), Comment("Extend the fit") };

    fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
  };

  class RegrowLoopHelix : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<RegrowLoopHelixConfig>;
      using KTRAJ = KinKal::LoopHelix;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PKTRAJPTR = std::unique_ptr<PKTRAJ>;
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
      using DOMAINPTR = std::shared_ptr<KinKal::Domain>;
      using DOMAINCOL = std::set<DOMAINPTR>;

      using KKMaterialConfig = KKMaterial::Config;

      explicit RegrowLoopHelix(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
      void endJob() override;
    private:
      int debug_;
      ProditionsHandle<StrawResponse> strawResponse_h_;
      ProditionsHandle<Tracker> alignedTracker_h_;
      std::unique_ptr<KinKal::BFieldMap> kkbf_;
      Config config_; // refit configuration object, containing the fit schedule
      KKFIT kkfit_;
      KKMaterial kkmat_;
      art::ProductToken<KalSeedCollection> kseedcol_T_;
      art::ProductToken<ComboHitCollection> chcol_T_;
      art::ProductToken<CaloClusterCollection> cccol_T_;
      art::ProductToken<IndexMap> indexmap_T_;
      art::InputTag ksmca_T_;
      bool fillMCAssns_;
      bool extend_;
      std::unique_ptr<KKExtrap> extrap_;
  };

  RegrowLoopHelix::RegrowLoopHelix(const Parameters& settings) : art::EDProducer(settings),
  debug_(settings().debug()),
  config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
  kkfit_(settings().kkfitSettings()),
  kkmat_(settings().matSettings()),
  kseedcol_T_(consumes<KalSeedCollection>(settings().kalSeedCollection())),
  chcol_T_(consumes<ComboHitCollection>(settings().comboHitCollection())),
  cccol_T_(mayConsume<CaloClusterCollection>(settings().caloClusterCollection())),
  indexmap_T_(consumes<IndexMap>(settings().indexMap())),
  fillMCAssns_(settings().kalSeedMCAssns(ksmca_T_)),
  extend_(settings().extend())
  {
    produces<KKTRKCOL>();
    produces<KalSeedCollection>();
    if(settings().extrapSettings())extrap_ = make_unique<KKExtrap>(*settings().extrapSettings(),kkmat_);
    if( fillMCAssns_){
      consumes<KalSeedMCAssns>(ksmca_T_);
      produces <KalSeedMCAssns>();
    }
  }

  void RegrowLoopHelix::beginRun(art::Run& run)
  {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void RegrowLoopHelix::produce(art::Event& event)
  {
    // proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    GeomHandle<mu2e::Tracker> nominalTracker_h;
    GeomHandle<Calorimeter> calo_h;
   // find input event data
    auto kseed_H = event.getValidHandle<KalSeedCollection>(kseedcol_T_);
    const auto& kseedcol = *kseed_H;
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    const auto& chcol = *ch_H;
    auto cc_H = event.getHandle<CaloClusterCollection>(cccol_T_);
    auto indexmap_H = event.getValidHandle<IndexMap>(indexmap_T_);
    const auto& indexmap = *indexmap_H;
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    art::Handle<KalSeedMCAssns> ksmca_H;
    // create outputs
    unique_ptr<KKTRKCOL> ktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> rgkseedcol(new KalSeedCollection );
    std::unique_ptr<KalSeedMCAssns> ksmca;
    // deal with MC
    if(fillMCAssns_){
      ksmca_H = event.getHandle<KalSeedMCAssns>(ksmca_T_);
      if(!ksmca_H)throw cet::exception("RECO")<<"mu2e::RegrowLoopHelix: No KalSeedMCAssns found" << endl;
      ksmca = std::unique_ptr<KalSeedMCAssns>(new KalSeedMCAssns);
    }
    size_t iseed(0);
    for (auto const& kseed : kseedcol) {
      if(!kseed.loopHelixFit())throw cet::exception("RECO")<<"mu2e::RegrowLoopHelix: passed KalSeed from non-LoopHelix fit " << endl;
      // regrow the components from the seed
      PKTRAJPTR trajptr = kseed.loopHelixFitTrajectory();
      KKSTRAWHITCOL strawhits;
      strawhits.reserve(kseed.hits().size());
      KKSTRAWXINGCOL strawxings;
      strawxings.reserve(kseed.straws().size());
      KKCALOHITCOL calohits;
      DOMAINCOL domains;
      // create the trajectory. This is done here to be strongly typed
      auto goodhits = kkfit_.regrowComponents(kseed, chcol, indexmap,
          *tracker,*calo_h,*strawresponse,*kkbf_, kkmat_.strawMaterial(),
          trajptr, strawhits, calohits, strawxings, domains);
      if(debug_ > 1){
        std::cout << "Regrew " << strawhits.size() << " straw hits, " << strawxings.size() << " straw xings, " << calohits.size() << " CaloHits and " << domains.size() << " domains, status = " << goodhits << std::endl;
      }
      if(debug_ > 2){
        unsigned nhactive(0);
        unsigned nsactive(0);
        for( auto const& strawh : strawhits)if(strawh->active())++nhactive;
        for( auto const& strawx : strawxings)if(strawx->active())++nsactive;
        std::cout << "Regrew " << nhactive << " active hits and " << nsactive << " active straws" << std::endl;
      }
      if(debug_ > 5)static_cast<KinKal::PiecewiseTrajectory<KTRAJ>*>(trajptr.get())->print(std::cout,2);
      // require hits and consistent BField domains
      if(goodhits && (domains.size() > 0 || !config_.bfcorr_)){
      // create the KKTrack from these components
        auto ktrk = std::make_unique<KKTRK>(config_,*kkbf_,kseed.particle(),trajptr,strawhits,strawxings,calohits,domains);
        if(ktrk && ktrk->fitStatus().usable()){
          if(debug_ > 0) std::cout << "RegrowLoopHelix: successful track refit" << std::endl;
          if(extend_)kkfit_.extendTrack(config_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H , *ktrk );
          if(ktrk->fitStatus().usable()){
            // extrapolate as requested
            if(extrap_)extrap_->extrapolate(*ktrk);
            // sample the fit as requested
            kkfit_.sampleFit(*ktrk);
            // convert to seed output format
            TrkFitFlag fitflag = kseed.status();
            fitflag.merge(TrkFitFlag::Regrown);
            auto rgks = kkfit_.createSeed(*ktrk,fitflag,*calo_h,*nominalTracker_h);
            rgkseedcol->push_back(rgks);
            if(fillMCAssns_){
              // find the MC assns
              auto ksmcai = (*ksmca_H)[iseed];
              auto origksp = art::Ptr<KalSeed>(kseed_H,iseed);
              // test this is the right ptr
              if(ksmcai.first != origksp)throw cet::exception("Reco")<<"mu2e::RegrowLoopHelix: wrong KalSeed ptr"<< std::endl;
              auto mcseedp = ksmcai.second;
              auto rgksp = art::Ptr<KalSeed>(KalSeedCollectionPID,rgkseedcol->size()-1,KalSeedCollectionGetter);
              ksmca->addSingle(rgksp,mcseedp);
              // add the original too
              ksmca->addSingle(origksp,mcseedp);
            }
            if(debug_ > 5)static_cast<const KinKal::PiecewiseTrajectory<KTRAJ>&>(ktrk->fitTraj()).print(std::cout,2);
            ktrkcol->push_back(ktrk.release());
          }
        } else if(debug_ > 0)
          std::cout << "RegrowLoopHelix: failed track refit, status " << ktrk->fitStatus() << std::endl;
      }
      ++iseed;
    }
    // store output
    event.put(move(ktrkcol));
    event.put(move(rgkseedcol));
    if(fillMCAssns_)event.put(move(ksmca));
  }

  void RegrowLoopHelix::endJob()
  {}

}
DEFINE_ART_MODULE(mu2e::RegrowLoopHelix)
