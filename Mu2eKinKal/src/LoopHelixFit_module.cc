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
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
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
#include "Offline/Mu2eKinKal/inc/ExtrapolateIPA.hh"
#include "Offline/Mu2eKinKal/inc/ExtrapolateST.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
// C++
#include <iostream>
#include <string>
#include <functional>
#include <vector>
#include <memory>
//
// Original author D. Brown (LBNL) 11/18/20
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
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KinKal::TimeDir;
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

  using KinKal::DVEC;
  using KinKal::VEC3;
  using MatEnv::DetMaterial;

  // extend the generic module configuration as needed
  struct KKLHModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections {Name("HelixSeedCollections"),     Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
  };
// Extrapolation configuration
  struct KKExtrapConfig {
    fhicl::Atom<float> Tol { Name("Tolerance"), Comment("Tolerance on fractional momemtum precision when extrapolating fits") };
    fhicl::Atom<float> MaxDt { Name("MaxDt"), Comment("Maximum time to extrapolate a fit") };
    fhicl::Atom<int> Debug { Name("Debug"), Comment("Debug level"), 0 };
  };
  struct LoopHelixFitConfig {
    fhicl::Table<KKLHModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
    fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
    fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
    fhicl::OptionalTable<KKFinalConfig> finalSettings { Name("FinalSettings") };
    fhicl::OptionalTable<KKExtrapConfig> Extrapolation { Name("Extrapolation") };
    // LoopHelix module specific config
    fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream") };
    fhicl::Atom<bool> pdgCharge { Name("UsePDGCharge"), Comment("Use particle charge from fitParticle")};
  };

  class LoopHelixFit : public art::EDProducer {
    public:
      using Parameters = art::EDProducer::Table<LoopHelixFitConfig>;
      explicit LoopHelixFit(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
    private:
      TrkFitFlag fitflag_;
      // parameter-specific functions that need to be overridden in subclasses
      KTRAJ makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const;
      bool goodFit(KKTRK const& ktrk) const;
      // extrapolation functions
      void extrapolate(KKTRK& ktrk) const;
      void extrapolateIPA(KKTRK& ktrk,TimeDir trkdir) const;
      void extrapolateST(KKTRK& ktrk,TimeDir trkdir) const;
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
      SurfaceMap smap_;
      bool extrapolate_;
      ExtrapolateToZ toTSDA_, toTracker_; // extrapolation predicate based on Z values
      ExtrapolateIPA extrapIPA_; // extrapolation to intersections with the IPA
      ExtrapolateST extrapST_; // extrapolation to intersections with the ST
      double ipathick_ = 0.511; // ipa thickness: should come from geometry service TODO
      double stthick_ = 0.1056; // st foil thickness: should come from geometry service TODO
  };

  LoopHelixFit::LoopHelixFit(const Parameters& settings) :  art::EDProducer{settings},
    fitflag_(TrkFitFlag::KKLoopHelix) ,
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
    fixedfield_(false),extrapolate_(false)
  {
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
    // configure extrapolation
    if(settings().Extrapolation()){
      extrapolate_ = true;
      auto const& IPA = smap_.DS().innerProtonAbsorberPtr();
      // global configs
      double maxdt = settings().Extrapolation()->MaxDt();
      double tol =  settings().Extrapolation()->Tol();
      int debug =  settings().Extrapolation()->Debug();
      // predicate to extrapolate through IPA
      extrapIPA_ = ExtrapolateIPA(maxdt,tol,IPA,debug);
      // predicate to extrapolate through ST
      std::cout << "IPA limits z " << extrapIPA_.zmin() << " " << extrapIPA_.zmax() << std::endl;
      extrapST_ = ExtrapolateST(maxdt,tol,smap_.ST(),debug);
      // temporary
      std::cout << "ST limits z " << extrapST_.zmin() << " " << extrapST_.zmax() << " r " << extrapST_.rmin() << " " << extrapST_.rmax() << std::endl;
      // extrapolate to the front of the tracker
      toTracker_ = ExtrapolateToZ(maxdt,tol,smap_.tracker().front().center().Z());
      // extrapolate to the back of the detector solenoid
      toTSDA_ = ExtrapolateToZ(maxdt,tol,smap_.DS().upstreamAbsorber().center().Z());
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
  }

  void LoopHelixFit::produce(art::Event& event ) {
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
            if(extrapolate_ && goodfit)extrapolate(*kktrk);
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

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const {
    auto const& helix = hseed.helix();
    DVEC pars;
    double psign = copysign(1.0,-charge*bnom.Z());
    pars[KTRAJ::rad_] = helix.radius()*psign;
    pars[KTRAJ::lam_] = helix.lambda()*fdir_.dzdt();
    pars[KTRAJ::cx_] = helix.centerx();
    pars[KTRAJ::cy_] = helix.centery();
    pars[KTRAJ::phi0_] = helix.fz0()+psign*M_PI_2;
    pars[KTRAJ::t0_] = hseed.t0().t0();
    // create the initial trajectory
    KinKal::Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    return KTRAJ(kkpars, mass_, charge, bnom, trange);
  }

  bool LoopHelixFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    bool retval = ktrk.fitStatus().usable() &&
      ktrk.fitTraj().front().parameterSign()*ktrk.seedTraj().front().parameterSign() > 0 &&
      ktrk.fitTraj().front().helicity()*ktrk.seedTraj().front().helicity() > 0;
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

  void LoopHelixFit::extrapolate(KKTRK& ktrk) const {
    // define the time direction according to the fit direction inside the tracker
    auto const& ftraj = ktrk.fitTraj();
    auto dir0 = ftraj.direction(ftraj.t0());
    TimeDir tdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    // extrapolate through the IPA in this direction
    extrapolateIPA(ktrk,tdir);
    // check if the particle has reflected
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto ipadir = ftraj.direction(endtime);
    if(ipadir.Z() * dir0.Z() > 0.0){
      // then extrapolate through the target
      extrapolateST(ktrk,tdir);
      // check if the particle has reflected
      double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
      auto stdir = ftraj.direction(endtime);
      if(stdir.Z() * dir0.Z() > 0.0){
        // extrapolate to the back of the DS
        ktrk.extrapolate(tdir,toTSDA_);
        // if we hit the TSDA we are done. Otherwise if we reflected, go back through the ST and IPA
        double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
        auto tsdadir = ftraj.direction(endtime);
        if(tsdadir.Z() * dir0.Z() < 0.0){
          extrapolateST(ktrk,tdir);
          extrapolateIPA(ktrk,tdir);
          ktrk.extrapolate(tdir,toTracker_);
        }
      } else {
        // reflection: extrapolate back through the IPA, then back to the tracker
        extrapolateIPA(ktrk,tdir);
        ktrk.extrapolate(tdir,toTracker_);
      }
    } else {
      // reflection: extrapolate back to the tracker entrance
      ktrk.extrapolate(tdir,toTracker_);
    }
  }

  void LoopHelixFit::extrapolateIPA(KKTRK& ktrk,TimeDir tdir) const {
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    extrapIPA_.reset();
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    do {
      ktrk.extrapolate(tdir,extrapIPA_);
      if(extrapIPA_.intersection().onsurface_ && extrapIPA_.intersection().inbounds_){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        auto const& IPA = smap_.DS().innerProtonAbsorberPtr();
        KKIPAXINGPTR ipaxingptr = std::make_shared<KKIPAXING>(IPA,IPASID,*kkmat_.IPAMaterial(),extrapIPA_.intersection(),reftrajptr,ipathick_,extrapIPA_.tolerance());
        if(extrapIPA_.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          ipaxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "IPA Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum() << std::endl;
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(extrapIPA_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapIPA_.intersection().onsurface_ && extrapIPA_.intersection().inbounds_);
  }

  void LoopHelixFit::extrapolateST(KKTRK& ktrk,TimeDir tdir) const {
    // extraplate the fit through the ST. This will add material effects for each foil intersection. It will continue till the
    // track exits the ST in Z
    extrapST_.reset();
    auto const& ftraj = ktrk.fitTraj();
    do {
      ktrk.extrapolate(tdir,extrapST_);
      if(extrapST_.intersection().onsurface_ && extrapST_.intersection().inbounds_){
        // we have a good intersection. Use this to create a Shell material Xing
        auto const& reftrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
        KKSTXINGPTR stxingptr = std::make_shared<KKSTXING>(extrapST_.foil(),extrapST_.foilId(),*kkmat_.STMaterial(),extrapST_.intersection(),reftrajptr,stthick_,extrapST_.tolerance());
        if(extrapST_.debug() > 0){
          double dmom, paramomvar, perpmomvar;
          stxingptr->materialEffects(dmom,paramomvar,perpmomvar);
          std::cout << "ST Xing dmom " << dmom << " para momsig " << sqrt(paramomvar) << " perp momsig " << sqrt(perpmomvar) << std::endl;
          std::cout << " before append mom = " << reftrajptr->momentum() << std::endl;
        }
        ktrk.addSTXing(stxingptr,tdir);
        if(extrapST_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapST_.intersection().onsurface_ && extrapST_.intersection().inbounds_);
  }

}
DEFINE_ART_MODULE(mu2e::LoopHelixFit)
