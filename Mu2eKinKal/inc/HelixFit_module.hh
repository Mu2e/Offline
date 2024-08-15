//
// KinKal fit module using an input helix seed.  Base for either central or loop helix fit
//
// Original author D. Brown (LBNL) 11/18/20
//
// framework
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
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/ExtraConfig.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
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
// root
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
namespace mu2e {
  using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
  using KKTRK = KKTrack<KTRAJ>;
  using KKTRKCOL = OwningPointerCollection<KKTRK>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;
  using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
  using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
  using KKSTRAWXING = KKStrawXing<KTRAJ>;
  using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
  using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
  using KKCALOHIT = KKCaloHit<KTRAJ>;
  using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
  using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KinKal::TimeDir;
  using KinKal::ExtraConfig;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eKinKal::KKFinalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;

  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;

  using KKMaterialConfig = KKMaterial::Config;
  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;

  // extend the generic module configuration as needed
  struct KKHelixModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections         {Name("HelixSeedCollections"),     Comment("Seed fit collections to be processed ") };
      fhicl::OptionalAtom<std::string> extrapDown { Name("ExtrapolateDownstreamSurface"), Comment("Extrapolate successful fits Downstream to this surface") };
      fhicl::OptionalAtom<std::string> extrapUp { Name("ExtrapolateUpstreamSurface"), Comment("Extrapolate successful fits Upstream to this surface") };
      fhicl::OptionalAtom<float> extrapTol { Name("ExtrapolationTolerance"), Comment("Tolerance on fractional momemtum precision when extrapolating fits") };
      fhicl::OptionalAtom<float> extrapMaxDt { Name("ExtrapolationMaxDt"), Comment("Maximum time to extrapolate a fit") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
  };

  struct HelixFitConfig {
    fhicl::Table<KKHelixModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::OptionalTable<KKFinalConfig> finalSettings { Name("FinalSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
    // helix module specific config
    fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream") };
    fhicl::Atom<bool> pdgCharge { Name("UsePDGCharge"), Comment("Use particle charge from fitParticle")};
  };

  class HelixFit : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<HelixFitConfig>;
      explicit HelixFit(const Parameters& settings,TrkFitFlag fitflag);
      virtual ~HelixFit() {}
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
    protected:
      TrkFitFlag fitflag_;
      // parameter-specific functions that need to be overridden in subclasses
      virtual KTRAJ makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const = 0;
      virtual bool goodFit(KKTRK const& ktrk) const = 0;
      // data payload
      std::vector<art::ProductToken<HelixSeedCollection>> hseedCols_;
      art::ProductToken<ComboHitCollection> chcol_T_;
      art::ProductToken<CaloClusterCollection> cccol_T_;
      TrkFitFlag goodseed_;
      bool saveall_;
      ProditionsHandle<StrawResponse> strawResponse_h_;
      ProditionsHandle<Tracker> alignedTracker_h_;
      int print_;
      PDGCode::type fpart_;
      TrkFitDirection fdir_;
      bool usePDGCharge_; // use the pdg particle charge: otherwise use the helicity and direction to determine the charge
      KKFIT kkfit_; // fit helper
      KKMaterial kkmat_; // material helper
      DMAT seedcov_; // seed covariance matrix
      double mass_; // particle mass
      int PDGcharge_; // PDG particle charge
      std::unique_ptr<KinKal::BFieldMap> kkbf_;
      Config config_; // initial fit configuration object
      Config exconfig_; // extension configuration object
      Config fconfig_; // final final configuration object
      bool fixedfield_; // special case usage for seed fits, if no BField corrections are needed
      ExtrapolateToZ extrap_; // extrapolation predicate based on Z values. Eventually this could be a more sophisticated test TODO
  };

  HelixFit::HelixFit(const Parameters& settings,TrkFitFlag fitflag) : art::EDProducer{settings},
    fitflag_(fitflag),
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    goodseed_(settings().modSettings().seedFlags()),
    saveall_(settings().modSettings().saveAll()),
    print_(settings().modSettings().printLevel()),
    fpart_(static_cast<PDGCode::type>(settings().modSettings().fitParticle())),
    fdir_(static_cast<TrkFitDirection::FitDirection>(settings().fitDirection())),
    usePDGCharge_(settings().pdgCharge()),
    kkfit_(settings().kkfitSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    fixedfield_(false) {
      // collection handling
      for(const auto& hseedtag : settings().modSettings().seedCollections()) { hseedCols_.emplace_back(consumes<HelixSeedCollection>(hseedtag)); }
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
      produces<KalHelixAssns>();
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::HelixFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      if(print_ > 0) std::cout << "Fit " << config_ << "Extension " << exconfig_;
      double bz(0.0);
      if(settings().modSettings().fixedBField(bz)){
        fixedfield_ = true;
        kkbf_ = std::move(std::make_unique<KKConstantBField>(VEC3(0.0,0.0,bz)));
      }
      // configure extrapolation
      std::string upsurf,downsur;
      if(settings().modSettings().extrapDown(downsurf) &&
          settings().modSettings().extrapUp(upsurf)){
        SurfaceMap smap;
        auto downsurf = std::dynamic_pointer_cast<KinKal::Plane>(smap.surface(downsurf));
        auto upsurf = std::dynamic_pointer_cast<KinKal::Plane>(smap.surface(upsurf));
        if((!downsurf) || (!upsurf) ||
            fabs(1.0 - downsurf->normal().Z()) > 1e-8 ||
            fabs(1.0 - upsurf->normal().Z()) > 1e-8 )
          throw cet::exception("RECO")<<"mu2e::HelixFit: invalid extrapolation surface ;must be planes orthogonal to z)" << endl;
        extrap_ = ExtrapolateToZ(settings().modSettings().extrapMaxDt(),
            settings().modSettings().extrapTol(),
            upsurf->center().Z(), downsurf->center().Z());
      }
      // setup optional fit finalization
      if(settings().finalSettings()){
        if(exconfig_.schedule_.size() > 0)
          fconfig_ = exconfig_;
        else
          fconfig_ = config_;
        fconfig_.schedule_.clear();
        fconfig_.schedule_.push_back(KinKal::MetaIterConfig()); // 1 stage, 0 temperature
        fconfig_.convdchisq_ = settings().finalSettings()->convdchisq();
        fconfig_.maxniter_ =  settings().finalSettings()->maxniter();
      }
    }

  void HelixFit::beginRun(art::Run& run) {
    // setup things that rely on data related to beginRun
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(fpart_).mass();
    PDGcharge_ = static_cast<int>(ptable->particle(fpart_).charge());
    // create KKBField
    if(!fixedfield_){
      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
    }
    if(print_ > 0) kkbf_->print(std::cout);
  }

  void HelixFit::produce(art::Event& event ) {
    // empty collections
    static MEASCOL nohits; // empty
    static EXINGCOL noexings; // empty
    // calo geom
    GeomHandle<Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getValidHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKTRKCOL> kktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
    unique_ptr<KalHelixAssns> kkseedassns(new KalHelixAssns());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // find the helix seed collections
    unsigned nseed(0);
    for (auto const& hseedtag : hseedCols_) {
      auto const& hseedcol_h = event.getValidHandle<HelixSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      nseed += hseedcol.size();
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];
        auto hptr = HPtr(hseedcol_h,iseed);
        // check helicity.  The test on the charge and helicity
        if(hseed.status().hasAllProperties(goodseed_) ){
          // test helix
          auto const& helix = hseed.helix();
          if(helix.radius() == 0.0 || helix.lambda() == 0.0 )
            throw cet::exception("RECO")<<"mu2e::HelixFit: degenerate seed parameters" << endl;
          auto zcent = Mu2eKinKal::zMid(hseed.hits());
          // take the magnetic field at the helix center as nominal
          VEC3 center(helix.centerx(), helix.centery(),zcent);
          auto bnom = kkbf_->fieldVect(center);
          // compute the charge from the helicity, fit direction, and BField direction
          double bz = bnom.Z();
          int charge = static_cast<int>(copysign(PDGcharge_,(-1)*helix.helicity().value()*fdir_.dzdt()*bz));
          // test consistency.  Modify this later when the HelixSeed knows which direction it's going TODO
          auto fitpart = fpart_;
          if(charge*PDGcharge_ < 0){
            if(usePDGCharge_)throw cet::exception("RECO")<<"mu2e::HelixFit: inconsistent charge" << endl;
            fitpart = static_cast<PDGCode::type>(-1*fitpart); // reverse sign
          }
          // time range of the hits
          auto trange = Mu2eKinKal::timeBounds(hseed.hits());
          // construt the seed trajectory
          KTRAJ seedtraj = makeSeedTraj(hseed,trange,bnom,charge);
          // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
          PKTRAJ pseedtraj(seedtraj);
          // first, we need to unwind the combohits.  We use this also to find the time range
          StrawHitIndexCollection strawHitIdxs;
          auto chcolptr = hseed.hits().fillStrawHitIndices(strawHitIdxs, StrawIdMask::uniquestraw);
          if(chcolptr != &chcol)
            throw cet::exception("RECO")<<"mu2e::KKHelixFit: inconsistent ComboHitCollection" << std::endl;
          // next, build straw hits and materials from these
          KKSTRAWHITCOL strawhits;
          strawhits.reserve(strawHitIdxs.size());
          KKSTRAWXINGCOL strawxings;
          strawxings.reserve(strawHitIdxs.size());
          if(kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings)){
            // optionally (and if present) add the CaloCluster as a constraint
            // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it TODO
            KKCALOHITCOL calohits;
            if (kkfit_.useCalo() && hseed.caloCluster().isNonnull())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
            // set the seed range given the hits and xings
            seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
            // create and fit the track
            auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits);
            // Check the fit
            auto goodfit = goodFit(*kktrk);
            // if we have an extension schedule, extend.
            if(goodfit && exconfig_.schedule().size() > 0) {
              kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *kktrk );
              goodfit = goodFit(*kktrk);
              // if finaling, apply that now.
              if(goodfit && fconfig_.schedule().size() > 0){
                kktrk->extend(fconfig_,nohits,noexings);
                goodfit = goodFit(*kktrk);
              }
            }
            // extrapolate as required
            if(goodfit && extrap_.size()>0) {
              // test the drection of this fit
              auto const& ftraj = kktrk->fitTraj();
              bool downstream = ftraj.momentum3(ftraj.range().mid()).Z() > 0.0; // replace with momentum at tracker middle TODO
              const static VEC3 opos(0.0,0.0,0.0);
              for(auto const& surf : extrap_){
                // configure the extrapolation time direction according to the surface and the track momentum direction
                if(surf.first.id() == SurfaceIdEnum::TT_Front){
                  xconfig_.xdir_ = downstream ? TimeDir::backwards : TimeDir::forwards;
                  double zpos = surf.second->tangentPlane(opos).center().Z(); // this is crude: I need an accessor that knows the TT_Front is a plane TODO
                  ExtrapolateToZ xtoz(*kktrk,xconfig_.xdir_,zpos);
                  kktrk->extrapolate(xconfig_,xtoz);
                } else if(surf.first.id() == SurfaceIdEnum::TT_Back){
                  xconfig_.xdir_ = downstream ? TimeDir::forwards : TimeDir::backwards;
                  double zpos = surf.second->tangentPlane(opos).center().Z();
                  ExtrapolateToZ xtoz(*kktrk,xconfig_.xdir_,zpos);
                  kktrk->extrapolate(xconfig_,xtoz);
                } else if(surf.first.id() == SurfaceIdEnum::TT_Outer){
                  // extrapolate in both time directions
                }
              }
            }

            if(print_>0)kktrk->printFit(std::cout,print_);
            if(goodfit || saveall_){
              TrkFitFlag fitflag(hptr->status());
              fitflag.merge(fitflag_);
              if(goodfit)
                fitflag.merge(TrkFitFlag::FitOK);
              else
                fitflag.clear(TrkFitFlag::FitOK);
              kkseedcol->push_back(kkfit_.createSeed(*kktrk,fitflag,*calo_h));
              // fill assns with the helix seed
              auto hptr = art::Ptr<HelixSeed>(hseedcol_h,iseed);
              auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
              kkseedassns->addSingle(kseedptr,hptr);
              // save (unpersistable) KKTrk in the event
              kktrkcol->push_back(kktrk.release());
            }
          }
        }
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << nseed << " Seeds" << std::endl;
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
  }

} // mu2e
