//
// KinKal fit module for LoopHelix
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
#include "Offline/Mu2eKinKal/inc/KKExtrap.hh"
// C++
#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
//
namespace mu2e {
  using KTRAJ= KinKal::LoopHelix;
  using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
  using KKTRK = KKTrack<KTRAJ>;
  using KKTRKCOL = OwningPointerCollection<KKTRK>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;
  using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
  using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
  using KKSTRAWXING = KKStrawXing<KTRAJ>;
  using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
  using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
  using KKSTXING = KKShellXing<KTRAJ,KinKal::Annulus>;
  using KKSTXINGPTR = std::shared_ptr<KKSTXING>;
  using KKSTXINGCOL = std::vector<KKSTXINGPTR>;
  using KKCALOHIT = KKCaloHit<KTRAJ>;
  using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
  using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KinKal::DVEC;
  using KinKal::TimeDir;
  using MatEnv::DetMaterial;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::Handle<CaloClusterCollection>;
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

  using CylPtr = std::shared_ptr<KinKal::Cylinder>;
  using DiskPtr = std::shared_ptr<KinKal::Disk>;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using FruPtr = std::shared_ptr<KinKal::Frustrum>;

  // extend the generic module configuration as needed
  struct KKLHModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections {Name("HelixSeedCollections"), Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
  };
  struct HelixMaskConfig {
    fhicl::OptionalAtom<float> minHelixP{ Name("MinHelixMom"), Comment("Minimum Momentum of a helix for a track to be fit.")};
  };
  struct LoopHelixFitConfig {
    fhicl::Table<KKLHModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
    fhicl::OptionalTable<KKFinalConfig> finalSettings { Name("FinalSettings") };
    fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
    // LoopHelix module specific config
    fhicl::OptionalAtom<double> slopeSigThreshold{ Name("SlopeSigThreshold"), Comment("Input helix seed slope significance threshold to assume the direction")};
    fhicl::OptionalAtom<std::string> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either \"upstream\" or \"downstream\"")};
    fhicl::Atom<bool> pdgCharge { Name("UsePDGCharge"), Comment("Use particle charge from fitParticle")};
    fhicl::OptionalTable<HelixMaskConfig> HelixMask { Name("HelixMask"), Comment("Selections applied to input helices")};
  };

  class LoopHelixFit : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<LoopHelixFitConfig>;
      explicit LoopHelixFit(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
      void endJob() override;
    private:
      TrkFitFlag fitflag_;
      // parameter-specific functions that need to be overridden in subclasses
      KTRAJ makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const;
      bool goodFit(KKTRK const& ktrk,KTRAJ const& seed) const;
      bool goodHelix(HelixSeed const& hseed) const;
      std::vector<TrkFitDirection::FitDirection> chooseHelixDir(HelixSeed const& hseed) const;
      std::unique_ptr<KKTRK> fitTrack(art::Event& event, HelixSeed const& hseed, const TrkFitDirection fdir, const PDGCode::type fitpart);
      void print_track_info(const KalSeed& kkseed, const KKTRK& ktrk) const;

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
      double slopeSigThreshold_; //helix slope significance threshold to assume the direction
      bool useHelixSlope_; //use the helix slope estimate to decide the fit direction
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
      std::unique_ptr<KKExtrap> extrap_; // extrapolation helper
      bool fixedfield_; // special case usage for seed fits, if no BField corrections are needed
      //Helix Mask params
      float minHelixP_ = -1.;
      int nSeen_ = 0;
      int nFit_ = 0;
      int nSkipped_ = 0;
      int nAmbiguous_ = 0;
      int nDownstream_ = 0;
      int nUpstream_ = 0;
  };

  LoopHelixFit::LoopHelixFit(const Parameters& settings) :  art::EDProducer{settings},
    fitflag_(TrkFitFlag::KKLoopHelix) ,
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    goodseed_(settings().modSettings().seedFlags()),
    saveall_(settings().modSettings().saveAll()),
    print_(settings().modSettings().printLevel()),
    fpart_(static_cast<PDGCode::type>(settings().modSettings().fitParticle())),
    useHelixSlope_(settings().slopeSigThreshold(slopeSigThreshold_)),
    usePDGCharge_(settings().pdgCharge()),
    kkfit_(settings().kkfitSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    fixedfield_(false)
    {
      std::string fdir;
      if(settings().fitDirection(fdir))fdir_ = fdir;
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
      // setup extrapolation
      if(settings().extrapSettings())extrap_ = make_unique<KKExtrap>(*settings().extrapSettings(),kkmat_);

      // setup optional fit finalization; this just updates the internals, not the fit result itself
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
      if (settings().HelixMask()){
        if (settings().HelixMask()->minHelixP())
          {minHelixP_ = settings().HelixMask()->minHelixP().value();}

      }
    }

  void LoopHelixFit::beginRun(art::Run& run) {
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

    // Print the fit direction configuration
    if(print_ > 0) printf("[LoopHelixFit::%s::%s] Fit dz/dt direction = %.0f, PDG = %i, use helix slope = %o with threshold %.1f\n",
        __func__, moduleDescription().moduleLabel().c_str(), fdir_.dzdt(), (int) fpart_,
        useHelixSlope_, (useHelixSlope_) ? slopeSigThreshold_ : -1.f);
  }

  std::vector<TrkFitDirection::FitDirection> LoopHelixFit::chooseHelixDir(HelixSeed const& hseed) const {
    std::vector<TrkFitDirection::FitDirection> fitdirs {};
    if(goodHelix(hseed)){
      if(fdir_ == TrkFitDirection::FitDirection::unknown)
        fitdirs = {TrkFitDirection::FitDirection::downstream, TrkFitDirection::FitDirection::upstream};
      else
        fitdirs = {fdir_.fitDirection()};
      // if using the helix slope to decide the fit direction, check its significance, and refine the list as needed
      if(useHelixSlope_) {
        auto predicted_dir = hseed.recoDir().predictDirection(slopeSigThreshold_);
        if(predicted_dir != TrkFitDirection::FitDirection::unknown){
          if(std::find(fitdirs.begin(), fitdirs.end(), predicted_dir) != fitdirs.end())
            fitdirs = {predicted_dir};
          else
            fitdirs.clear();
        }
      }
    }
    return fitdirs;
  }

  std::unique_ptr<KKTRK> LoopHelixFit::fitTrack(art::Event& event, HelixSeed const& hseed, const TrkFitDirection fdir, PDGCode::type fitpart) {
    // check the input
    if(fdir.fitDirection() != TrkFitDirection::FitDirection::downstream && fdir.fitDirection() != TrkFitDirection::FitDirection::upstream)
      throw cet::exception("RECO") << "mu2e::LoopHelixFit: Unknown helix propagation direction " << fdir.name();

    // Retrieve event information
    // calo geom
    GeomHandle<Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // empty collections
    static MEASCOL nohits; // empty
    static EXINGCOL noexings; // empty

    // Initialize helix-specific information
    auto const& helix = hseed.helix();
    auto zcent = Mu2eKinKal::zMid(hseed.hits());
    // take the magnetic field at the helix center as nominal
    VEC3 center(helix.centerx(), helix.centery(),zcent);
    static const double rhomax = 700.0; // this should come from conditions
    if(center.Rho() > rhomax) center = VEC3(rhomax*cos(center.Phi()),rhomax*sin(center.Phi()),center.Z());
    auto bnom = kkbf_->fieldVect(center);
    // compute the charge from the helicity, fit direction, and BField direction
    double bz = bnom.Z();
    const int charge = static_cast<int>(copysign(PDGcharge_,(-1)*helix.helicity().value()*bz*fdir.dzdt()));
    // test consistency.  Modify this later when the HelixSeed knows which direction it's going TODO
    if(charge*PDGcharge_ < 0){
      if(usePDGCharge_)throw cet::exception("RECO")<<"mu2e::HelixFit: inconsistent charge" << endl;
      fitpart = static_cast<PDGCode::type>(-1*fitpart); // reverse sign
    }
    // time range of the hits
    auto trange = Mu2eKinKal::timeBounds(hseed.hits());
    // Begin constructing the track fit
    // construct the seed trajectory
    KTRAJ seedtraj = makeSeedTraj(hseed,trange,bnom,charge);
    // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
    PTRAJ pseedtraj(seedtraj);
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
    if(!kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings)) {
      if(print_>0) printf("[LoopHelixFit::%s] Failed to create a track\n", __func__);
      return nullptr;
    }
    // optionally (and if present) add the CaloCluster as a constraint
    // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it TODO
    KKCALOHITCOL calohits;
    if (kkfit_.useCalo() && hseed.caloCluster().isNonnull()) {
      kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
    }
    // set the seed range given the hits and xings
    seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
    // create and fit the track
    auto ktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits);
    if(!ktrk) // check that the track exists
      throw cet::exception("RECO")<<"mu2e::LoopHelixFit: Track fit was performed but no track is found\n";

    auto goodfit = goodFit(*ktrk,seedtraj);
    if(print_>0) printf("[LoopHelixFit::%s] Before extending the fit: goodFit = %o, fitcon = %.4f, nHits = %2lu, %lu calo-hits\n",
        __func__, goodfit, ktrk->fitStatus().chisq_.probability(), ktrk->strawHits().size(), ktrk->caloHits().size());
    // if we have an extension schedule, extend.
    if(goodfit && exconfig_.schedule().size() > 0) {
      kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *ktrk );
      goodfit = goodFit(*ktrk,seedtraj);
      // if finaling, apply that now.
      if(goodfit && fconfig_.schedule().size() > 0){
        ktrk->extend(fconfig_,nohits,noexings);
        goodfit = goodFit(*ktrk,seedtraj);
      }
    }
    if(print_>0) printf("[LoopHelixFit::%s] After extending the fit : goodFit = %o, fitcon = %.4f, nHits = %2lu, %lu calo-hits\n",
        __func__, goodfit, ktrk->fitStatus().chisq_.probability(), ktrk->strawHits().size(), ktrk->caloHits().size());
    if((!goodfit) && (! saveall_)) ktrk.reset();
    return ktrk;
  }

  void LoopHelixFit::produce(art::Event& event ) {
    // calo geom
    GeomHandle<Calorimeter> calo_h;
    // create output
    unique_ptr<KKTRKCOL> ktrkcol(new KKTRKCOL );
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
        ++nSeen_;
        auto const& hseed = hseedcol[iseed];
        // determine the fit direction hypotheses + check whether helix momentum is over the minimum momentum threshold (if set)
        auto helix_dirs = chooseHelixDir(hseed);
        if(helix_dirs.empty()) {
          ++nSkipped_;
          continue; //bad helix, no fits to perform
        }
        ++nFit_;
        const unsigned dirs_size = helix_dirs.size();
        const bool undefined_dir = dirs_size > 1; //fitting multiple hypotheses to determine the best fit
        if(undefined_dir) ++nAmbiguous_;
        // fit each track hypothesis
        for(auto helix_dir : helix_dirs) {
          auto ktrk = fitTrack(event, hseed,  TrkFitDirection(helix_dir), fpart_);
          if(!ktrk) continue; //ensure that the track exists
          // extrapolate as required
          if(extrap_)extrap_->extrapolate(*ktrk);
          if(print_>1) ktrk->printFit(std::cout,print_-1);
          auto hptr = HPtr(hseedcol_h,iseed);
          TrkFitFlag fitflag(hptr->status());
          fitflag.merge(fitflag_);
          if(undefined_dir) fitflag.merge(TrkFitFlag::AmbFitDir);
          // sample the fit as requested
          kkfit_.sampleFit(*ktrk);
          // convert to seed output format
          auto kkseed = kkfit_.createSeed(*ktrk,fitflag,*calo_h);
          if(print_>0) print_track_info(kkseed, *ktrk);
          kkseedcol->push_back(kkseed);
          // fill assns with the helix seed
          auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
          kkseedassns->addSingle(kseedptr,hptr);
          // save (unpersistable) KKTrk in the event
          ktrkcol->push_back(ktrk.release());
          //increment the counts
          if(helix_dir == TrkFitDirection::FitDirection::downstream) ++nDownstream_;
          if(helix_dir == TrkFitDirection::FitDirection::upstream  ) ++nUpstream_;
        } //end track fit result loop
      } //end helix seed loop
    } //end helix colllection loop

    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << ktrkcol->size() << " tracks from " << nseed << " Seeds" << std::endl;
    event.put(move(ktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
  }

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const {
    auto const& helix = hseed.helix();
    DVEC pars;
    double psign = copysign(1.0,-charge*bnom.Z());
    pars[KTRAJ::rad_] = helix.radius()*psign;
    pars[KTRAJ::lam_] = helix.lambda();
    pars[KTRAJ::cx_] = helix.centerx();
    pars[KTRAJ::cy_] = helix.centery();
    pars[KTRAJ::phi0_] = helix.fz0()+psign*M_PI_2;
    pars[KTRAJ::t0_] = hseed.t0().t0();
    // create the initial trajectory
    KinKal::Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    KTRAJ ktraj(kkpars, mass_, charge, bnom, trange);
    // test position and direction and z=0
    auto hdir = helix.direction(0.0);
    auto kdir = ktraj.direction(ktraj.t0());
    double dirdot = hdir.Dot(kdir);
    // the original helix doesn't have a time direction (geometric helix) so allow both interpretations
    // the tolerance in the test allows for a difference between global and local parameters (B not along Z axis)
    if(1.0- fabs(dirdot) > 1e-2)throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Seed helix translation error, dirdot = " << dirdot << endl;
    return ktraj;
  }

  bool LoopHelixFit::goodFit(KKTRK const& ktrk,KTRAJ const& seed) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity. Test at the t0 segment
    auto t0 = Mu2eKinKal::zTime(ktrk.fitTraj(),0.0,ktrk.fitTraj().range().mid());
    auto const& t0seg = ktrk.fitTraj().nearestPiece(t0);
    bool retval = ktrk.fitStatus().usable() && t0seg.parameterSign()*seed.parameterSign() > 0 && t0seg.helicity()*seed.helicity() > 0;
    // also check that the fit is inside the physical detector volume.  Test where the StrawHits are
    if(retval){
      for(auto const& shptr : ktrk.strawHits()) {
        if(shptr->active() && !Mu2eKinKal::inDetector(ktrk.fitTraj().position3(shptr->time()))){
          retval = false;
          break;
        }
      }
    }
    // test that the trajectory is inside the DS
    if(retval){
      static unsigned ntimes(100);
      double dt = ktrk.fitTraj().range().range()/(ntimes-1);
      for(unsigned it=0;it< ntimes; ++it) {
        double ttest = ktrk.fitTraj().range().begin() + it*dt;
        auto tpos = ktrk.fitTraj().position3(ttest);
        if(!ktrk.bfield().inRange(tpos)){
          retval = false;
          break;
        }
      }
    }
    return retval;
  }

  bool LoopHelixFit::goodHelix(HelixSeed const& hseed) const {
    // check helicity.  The test on the charge and helicity
    if(!hseed.status().hasAllProperties(goodseed_)) return false;
    // test helix
    auto const& helix = hseed.helix();
    if(helix.radius() == 0.0 || helix.lambda() == 0.0) {
      if(print_>0) printf("LoopHelixFit::%s::%s: Degenerate helix seed parameters, skipping helix\n", __func__, moduleDescription().moduleLabel().c_str());
      return false;
    }
    auto zcent = Mu2eKinKal::zMid(hseed.hits());
    // take the magnetic field at the helix center as nominal
    VEC3 center(helix.centerx(), helix.centery(),zcent);
    auto bnom = kkbf_->fieldVect(center);
    // convert momentum from mm to MeV
    const double bz = bnom.Z();
    const double b_to_p(3./10.); // conversion factor
    if (std::fabs(helix.momentum()*bz*b_to_p) < minHelixP_) {
      if(print_>0) printf("LoopHelixFit::%s::%s: Helix does not pass minimum momentum threshold. \n", __func__, moduleDescription().moduleLabel().c_str());
      return false;
    }
    return true;
  }

  void LoopHelixFit::print_track_info(const KalSeed& kkseed, const KKTRK& ktrk) const {
    if(print_>0) printf("[LoopHelixFit::%s::%s] Create seed             : fitcon = %.4f, nHits = %2lu, seedActiveHits = %2u, %lu calo-hits\n",
        __func__, moduleDescription().moduleLabel().c_str(), ktrk.fitStatus().chisq_.probability(), ktrk.strawHits().size(),
        kkseed.nHits(), ktrk.caloHits().size());
    if(print_>1) { //print the hit flags for the track and the KalSeed as well as the KalSegments and intersections
      for(auto const& hit : ktrk.strawHits())
        printf("  [LoopHelixFit::%s::%s] KKTRK straw flags: %s", __func__, moduleDescription().moduleLabel().c_str(), hit->hit().flag().hex().c_str());
      for(auto const& hit : kkseed.hits())
        printf("  [LoopHelixFit::%s::%s] Seed straw flags : %s", __func__, moduleDescription().moduleLabel().c_str(), hit.flag().hex().c_str());
      printf("[LoopHelixFit::%s::%s] KalSeed has %lu intersections\n", __func__, moduleDescription().moduleLabel().c_str(), kkseed.intersections().size());
      for(size_t ikinter = 0; ikinter < kkseed.intersections().size(); ++ikinter){
        auto const& kinter = kkseed.intersections()[ikinter];
        printf("[LoopHelixFit::%s::%s]   Seed %10s intersection: t0 = %6.1f, t0Err = %6.4f, mom = %6.2f, momErr = %6.4f\n",
            __func__, moduleDescription().moduleLabel().c_str(), kinter.surfaceId().name().c_str(), kinter.time(), std::sqrt(kinter.loopHelix().paramVar(KinKal::LoopHelix::t0_)),
            kinter.mom(), kinter.momerr());
      }
      printf("[LoopHelixFit::%s::%s] KalSeed has %lu segments\n", __func__, moduleDescription().moduleLabel().c_str(), kkseed.segments().size());
      for(size_t ikseg = 0; ikseg < kkseed.segments().size(); ++ikseg){
        auto const& kseg = kkseed.segments()[ikseg];
        printf("[LoopHelixFit::%s::%s]   Seed segment %lu: tmin = %6.1f, terr = %6.4f, mom = %6.2f, momErr = %6.4f\n",
            __func__, moduleDescription().moduleLabel().c_str(), ikseg, kseg.tmin(), std::sqrt(kseg.loopHelix().paramVar(KinKal::LoopHelix::t0_)),
            kseg.mom(), kseg.momerr());
      }
    }

  }

  void LoopHelixFit::endJob() {
    if(print_ > 0) {
      printf("[LoopHelixFit::%s::%s] Saw %i helix seeds, %i had ambiguous dz/dt slopes, accepted %i downstream and %i upstream fits\n",
                          __func__, moduleDescription().moduleLabel().c_str(), nSeen_, nAmbiguous_, nDownstream_, nUpstream_);
      printf("Number of fits: %i;  number of helices skipped: %i \n ", nFit_, nSkipped_);
    }
  }
}

DEFINE_ART_MODULE(mu2e::LoopHelixFit)
