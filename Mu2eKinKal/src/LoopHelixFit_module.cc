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
  using KinKal::DVEC;
  using KinKal::TimeDir;
  using MatEnv::DetMaterial;
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

  using DiskPtr = std::shared_ptr<KinKal::Disk>;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using FruPtr = std::shared_ptr<KinKal::Frustrum>;

  // extend the generic module configuration as needed
  struct KKLHModuleConfig : KKModuleConfig {
    fhicl::Sequence<art::InputTag> seedCollections {Name("HelixSeedCollections"),     Comment("Seed fit collections to be processed ") };
    fhicl::OptionalAtom<double> fixedBField { Name("ConstantBField"), Comment("Constant BField value") };
    fhicl::Sequence<std::string> sampleSurfaces { Name("SampleSurfaces"), Comment("When creating the KalSeed, sample the fit at these surfaces") };
    fhicl::Atom<bool> sampleInRange { Name("SampleInRange"), Comment("Require sample times to be inside the fit trajectory time range") };
    fhicl::Atom<bool> sampleInBounds { Name("SampleInBounds"), Comment("Require sample intersection point be inside surface bounds (within tolerance)") };
    fhicl::Atom<float> sampleTol { Name("SampleTolerance"), Comment("Tolerance for sample surface intersections (mm)") };
  };
  // Extrapolation configuration
  struct KKExtrapConfig {
    fhicl::Atom<float> Tol { Name("Tolerance"), Comment("Tolerance on fractional momemtum precision when extrapolating fits") };
    fhicl::Atom<float> MaxDt { Name("MaxDt"), Comment("Maximum time to extrapolate a fit") };
    fhicl::Atom<bool> BackToTracker { Name("BackToTracker"), Comment("Extrapolate reflecting tracks back to the tracker") };
    fhicl::Atom<bool> ToTrackerEnds { Name("ToTrackerEnds"), Comment("Extrapolate tracks to the tracker ends") };
    fhicl::Atom<bool> Upstream { Name("Upstream"), Comment("Extrapolate tracks upstream") };
    fhicl::Atom<bool> ToOPA { Name("ToOPA"), Comment("Test tracks for intersection with the OPA") };
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
    // fhicl::Atom<bool> useHelixSlope{ Name("UseHelixSlope"), Comment("Use the helix slope to decide the particle fit direction (upstream or downstream)"), false };
    // fhicl::Atom<double> slopeSigThreshold{ Name("SlopeSigThreshold"), Comment("Helix slope significance threshold to assume the direction"), [this](){ return useHelixSlope(); }};
    fhicl::OptionalAtom<double> slopeSigThreshold{ Name("SlopeSigThreshold"), Comment("Helix slope significance threshold to assume the direction")};
    fhicl::OptionalAtom<double> slopeSigCut{ Name("SlopeSigCut"), Comment("Helix slope significance cut when assuming a fit direction")};
    fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream"), [this](){ return !slopeSigThreshold(); }};
    fhicl::Atom<bool> pdgCharge { Name("UsePDGCharge"), Comment("Use particle charge from fitParticle")};
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
      bool goodFit(KKTRK const& ktrk) const;
      // extrapolation functions
      void extrapolate(KKTRK& ktrk) const;
      void toTrackerEnds(KKTRK& ktrk,VEC3 const& dir0) const;
      bool extrapolateIPA(KKTRK& ktrk,TimeDir trkdir) const;
      bool extrapolateST(KKTRK& ktrk,TimeDir trkdir) const;
      bool extrapolateTracker(KKTRK& ktrk,TimeDir tdir) const;
      bool extrapolateTSDA(KKTRK& ktrk,TimeDir tdir) const;
      void toOPA(KKTRK& ktrk, double tstart, TimeDir tdir) const;
      void sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const;

      // FIXME: Remove debug function
      unsigned activeHits(KKTRK& ktrk) const {
        unsigned nactive(0);
        for(auto const& hit : ktrk.strawHits()) {
          if(hit->hit().flag().hasAllProperties(StrawHitFlag::active)) ++nactive;
        }
        return nactive;
      }
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
      double slopeSigCut_; //helix slope significance to cut on when assuming a fit direction
      bool useSlopeSigCut_; //apply a helix slope significance cut
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
      double sampletol_; // surface intersection tolerance (mm)
      bool sampleinrange_, sampleinbounds_; // require samples to be in range or on surface
      bool fixedfield_; // special case usage for seed fits, if no BField corrections are needed
      SurfaceMap smap_;
      AnnPtr tsdaptr_;
      DiskPtr trkfrontptr_, trkmidptr_, trkbackptr_;
      FruPtr opaptr_;
      bool extrapolate_, backToTracker_, toOPA_, toTrackerEnds_, upstream_;
      ExtrapolateToZ TSDA_, trackerFront_, trackerBack_; // extrapolation predicate based on Z values
      ExtrapolateIPA extrapIPA_; // extrapolation to intersections with the IPA
      ExtrapolateST extrapST_; // extrapolation to intersections with the ST
      double ipathick_ = 0.511; // ipa thickness: should come from geometry service TODO
      double stthick_ = 0.1056; // st foil thickness: should come from geometry service TODO
      SurfaceMap::SurfacePairCollection sample_; // surfaces to sample the fit
      int nSeen_ = 0;
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
    useSlopeSigCut_(settings().slopeSigCut(slopeSigCut_)),
    fdir_(static_cast<TrkFitDirection::FitDirection>((useHelixSlope_) ? 0 : settings().fitDirection())),
    usePDGCharge_(settings().pdgCharge()),
    kkfit_(settings().kkfitSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    sampletol_(settings().modSettings().sampleTol()),
    sampleinrange_(settings().modSettings().sampleInRange()),
    sampleinbounds_(settings().modSettings().sampleInBounds()),
    fixedfield_(false), extrapolate_(false), backToTracker_(false), toOPA_(false)
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
      if(useSlopeSigCut_ && useHelixSlope_) //can only cut on significance if using a fixed fit direction hypothesis
        throw cet::exception("RECO")<<"mu2e::LoopHelixFit: Configuration error. Using helix slope for fit direction and cutting on the slope!"<< endl;
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
        backToTracker_ = settings().Extrapolation()->BackToTracker();
        toTrackerEnds_ = settings().Extrapolation()->ToTrackerEnds();
        upstream_ = settings().Extrapolation()->Upstream();
        toOPA_ = settings().Extrapolation()->ToOPA();
        auto const& IPA = smap_.DS().innerProtonAbsorberPtr();
        // global configs
        double maxdt = settings().Extrapolation()->MaxDt();
        double tol =  settings().Extrapolation()->Tol();
        int debug =  settings().Extrapolation()->Debug();
        // predicate to extrapolate through IPA
        extrapIPA_ = ExtrapolateIPA(maxdt,tol,IPA,debug);
        // predicate to extrapolate through ST
        if(debug > 0)std::cout << "IPA limits z " << extrapIPA_.zmin() << " " << extrapIPA_.zmax() << std::endl;
        extrapST_ = ExtrapolateST(maxdt,tol,smap_.ST(),debug);
        // temporary
        if(debug > 0)std::cout << "ST limits z " << extrapST_.zmin() << " " << extrapST_.zmax() << " r " << extrapST_.rmin() << " " << extrapST_.rmax() << std::endl;
        // extrapolate to the front of the tracker
        trackerFront_ = ExtrapolateToZ(maxdt,tol,smap_.tracker().front().center().Z(),debug);
        trackerBack_ = ExtrapolateToZ(maxdt,tol,smap_.tracker().back().center().Z(),debug);
        // extrapolate to the back of the detector solenoid
        TSDA_ = ExtrapolateToZ(maxdt,tol,smap_.DS().upstreamAbsorber().center().Z(),debug);
        tsdaptr_ = smap_.DS().upstreamAbsorberPtr();
        trkfrontptr_ = smap_.tracker().frontPtr();
        trkmidptr_ = smap_.tracker().middlePtr();
        trkbackptr_ = smap_.tracker().backPtr();
        opaptr_ = smap_.DS().outerProtonAbsorberPtr();
      }

      // additional surfaces to sample: these should be replaced by extrapolation TODO
      SurfaceIdCollection ssids;
      for(auto const& sidname : settings().modSettings().sampleSurfaces()){
        ssids.push_back(SurfaceId(sidname,-1)); // match all elements
      }
      // translate the sample and extend surface names to actual surfaces using the SurfaceMap.  This should come from the
      // geometry service eventually, TODO
      SurfaceMap smap;
      smap.surfaces(ssids,sample_);
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
    printf("[LoopHelixFit::%s::%s] Fit dz/dt direction = %.0f, PDG = %i, use helix slope = %o with threshold %.1f, use helix slope cut = %o with theshold %.1f\n",
           __func__, moduleDescription().moduleLabel().c_str(), fdir_.dzdt(), (int) fpart_,
           useHelixSlope_, (useHelixSlope_) ? slopeSigThreshold_ : -1.f,
           useSlopeSigCut_, (useSlopeSigCut_) ? slopeSigCut_ : -999.f);
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
        auto hptr = HPtr(hseedcol_h,iseed);
        // check helicity.  The test on the charge and helicity
        if(!hseed.status().hasAllProperties(goodseed_)) continue;
        // test helix
        auto const& helix = hseed.helix();
        if(helix.radius() == 0.0 || helix.lambda() == 0.0)
          throw cet::exception("RECO")<<"mu2e::HelixFit: degenerate seed parameters" << endl;
        // test the helix slope significance if requested
        if(useSlopeSigCut_) {
          const float slope    = hseed.recoDir().slope();
          const float slopeErr = std::fabs(hseed.recoDir().slopeErr());
          const float slopeSig = (slopeErr > 0.f) ? slope/slopeErr : 0.f;
          if(slopeSig*fdir_.dzdt() < slopeSigCut_) { //slope is below the given threshold relative to the given direction
            if(print_ > 0) printf("[LoopHelixFit::%s] Skipping helix seed with slope significance %.1f\n", __func__, slopeSig);
            continue;
          }
        }
        auto zcent = Mu2eKinKal::zMid(hseed.hits());
        // take the magnetic field at the helix center as nominal
        VEC3 center(helix.centerx(), helix.centery(),zcent);
        static const double rhomax = 700.0; // this should come from conditions
        if(center.Rho() > rhomax) center = VEC3(rhomax*cos(center.Phi()),rhomax*sin(center.Phi()),center.Z());
        auto bnom = kkbf_->fieldVect(center);
        // compute the charge from the helicity, fit direction, and BField direction
        double bz = bnom.Z();
        const int helix_dir = (useHelixSlope_) ? hseed.recoDir().predictDirection(slopeSigThreshold_) : int(fdir_.dzdt());
        const bool undefined_dir = helix_dir == HelixRecoDir::PropDir::ambiguous; //check if no direction is significantly preferred
        if(undefined_dir) ++nAmbiguous_;
        double pchi2_seed [] = {-1., -1.}; //p(chi^2,N(dof)) for {downstream, upstream} hypothesis at the seed fit stage
        double pchi2_drift[] = {-1., -1.}; // at the drift fit stage
        std::vector<std::unique_ptr<KKTRK>> ktrks;
        const int nloops = (undefined_dir) ? 2 : 1;
        for(int iloop = 0; iloop < nloops; ++iloop) {  //if an undefined direction, check both directions, pick the best fit result
          const float dir = (undefined_dir) ? -2.f*iloop + 1.f : 1.f*helix_dir;
          const int charge = static_cast<int>(copysign(PDGcharge_,(-1)*helix.helicity().value()*dir*bz));
          // test consistency.  Modify this later when the HelixSeed knows which direction it's going TODO
          auto fitpart = fpart_;
          if(charge*PDGcharge_ < 0){
            if(usePDGCharge_)throw cet::exception("RECO")<<"mu2e::HelixFit: inconsistent charge" << endl;
            fitpart = static_cast<PDGCode::type>(-1*fitpart); // reverse sign
          }
          // time range of the hits
          auto trange = Mu2eKinKal::timeBounds(hseed.hits());
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
          if(kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings)){
            // optionally (and if present) add the CaloCluster as a constraint
            // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it TODO
            KKCALOHITCOL calohits;
            if (kkfit_.useCalo() && hseed.caloCluster().isNonnull()) {
              kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
            }
            // set the seed range given the hits and xings
            seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
            // create and fit the track
            if(print_>0 && undefined_dir) printf("[LoopHelixFit::%s] Creating track fit direction hypothesis %i\n", __func__, iloop);
            ktrks.push_back(make_unique<KKTRK>(config_,*kkbf_,seedtraj,fitpart,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits));
            auto& ktrk = ktrks[iloop];
            auto goodfit = goodFit(*ktrk);
            if(goodfit) pchi2_seed[iloop] = ktrk->fitStatus().chisq_.probability(); //store the chi^2 before extending the fit
            if(print_>0) printf("[LoopHelixFit::%s] Before extending the fit: fitcon = %.4f, nHits = %2lu, %lu calo-hits\n",
                                __func__, pchi2_seed[iloop], ktrk->strawHits().size(), ktrk->caloHits().size());
            // if we have an extension schedule, extend.
            if(goodfit && exconfig_.schedule().size() > 0) {
              kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *ktrk );
              goodfit = goodFit(*ktrk);
              // if finaling, apply that now.
              if(goodfit && fconfig_.schedule().size() > 0){
                ktrk->extend(fconfig_,nohits,noexings);
                goodfit = goodFit(*ktrk);
              }
            }

            //store the fit quality result if it's a good fit
            if(goodfit) pchi2_drift[iloop] = ktrk->fitStatus().chisq_.probability();
            if(print_>0) printf("[LoopHelixFit::%s] After extending the fit : fitcon = %.4f, nHits = %2lu, %lu calo-hits\n",
                                __func__, pchi2_drift[iloop], ktrk->strawHits().size(), ktrk->caloHits().size());
            if(print_>0 && undefined_dir)
              printf("[LoopHelixFit::%s] Track fit hypothesis %i results: Status = %i, p(chi^2) = %.4f (%.4f at seed fit), N(hits) = %2lu, N(calo hits) = %lu, goodFit = %o\n",
                     __func__, iloop, ktrk->fitStatus().status_, pchi2_drift[iloop], pchi2_seed[iloop],
                     ktrk->strawHits().size(), ktrk->caloHits().size(), goodfit);
          } else {
            ktrks.push_back(nullptr);
          }
        }
        //determine the best fit if needed
        int index = 0;
        if(undefined_dir) {
          //if only one fit result is undefined, use the other fit (if both are undefined it doesn't matter which we pick)
          if(!ktrks[0] || pchi2_drift[0] < 0.) index = 1;
          else if(!ktrks[1] || pchi2_drift[1] < 0.) index = 0;
          //check if a calo-cluster was dropped for only one of the fits, use the seed fit instead as dropping the cluster is a bad sign
          else if(ktrks[0]->caloHits().size() != ktrks[1]->caloHits().size()) {
            //pick the fit that included the cluster in the result
            index = (ktrks[0]->caloHits().size() > ktrks[1]->caloHits().size()) ? 0 : 1;
          } else index = (pchi2_drift[0] >= pchi2_drift[1]) ? 0 : 1; //default to using the full fit p(chi^2) results
        }

        if(pchi2_drift[index] < 0.) continue; //failed earlier tests

        if(print_ > 0 && undefined_dir)
          printf("[LoopHelixFit::%s]: Fit both direction options, slope = %9.2e, sig = %.2f, status down = %i, status up = %i, p(chi^2) down = %.4f, p(chi^2) up = %.4f: chose %s\n",
                 __func__, hseed.recoDir().slope(), hseed.recoDir().slopeSig(),
                 ktrks[0] && goodFit(*ktrks[0]), ktrks[1] && goodFit(*ktrks[1]),
                 pchi2_drift[0], pchi2_drift[1], (index == 0) ? "downstream" : "upstream");

        if(!ktrks[index]) //ensure that the track was properly created
          throw cet::exception("RECO")<<"mu2e::LoopHelixFit: Track fit was performed but no track is found\n";

        auto& ktrk = ktrks[index];

        // Check the fit
        auto goodfit = goodFit(*ktrk);

        // extrapolate as required
        if(goodfit && extrapolate_) extrapolate(*ktrk);
        if(print_>1) ktrk->printFit(std::cout,print_-1);
        if(goodfit || saveall_){
          TrkFitFlag fitflag(hptr->status());
          fitflag.merge(fitflag_);
          if(goodfit)
            fitflag.merge(TrkFitFlag::FitOK);
          else
            fitflag.clear(TrkFitFlag::FitOK);
          if(undefined_dir) fitflag.merge(TrkFitFlag::AmbFitDir);
          auto kkseed = kkfit_.createSeed(*ktrk,fitflag,*calo_h);
          if(print_>0) printf("[LoopHelixFit::%s] Create seed             : fitcon = %.4f, nHits = %2lu, seedActiveHits = %2u, %lu calo-hits\n",
                              __func__, ktrk->fitStatus().chisq_.probability(), ktrk->strawHits().size(),
                              kkseed.nHits(), ktrk->caloHits().size());
          if(print_>1) { //print the hit flags for the track and the KalSeed
            for(auto const& hit : ktrk->strawHits())
              printf("  [LoopHelixFit::%s] KKTRK straw flags: %s", __func__, hit->hit().flag().hex().c_str());
            for(auto const& hit : kkseed.hits())
              printf("  [LoopHelixFit::%s] Seed straw flags : %s", __func__, hit.flag().hex().c_str());
          }
          sampleFit(*ktrk,kkseed._inters);
          if(print_>0) {
            printf("[LoopHelixFit::%s] KalSeed has %lu intersections\n", __func__, kkseed.intersections().size());
            for(size_t ikinter = 0; ikinter < kkseed.intersections().size(); ++ikinter){
              auto const& kinter = kkseed.intersections()[ikinter];
              printf("[LoopHelixFit::%s] Seed %10s intersection: t0 = %.1f, t0Err = %.4f, mom = %.2f, momErr = %.4f\n",
                     __func__, kinter.surfaceId().name().c_str(), kinter.time(), std::sqrt(kinter.loopHelix().paramVar(KinKal::LoopHelix::t0_)),
                     kinter.mom(), kinter.momerr());
            }
          }
          kkseedcol->push_back(kkseed);
          // fill assns with the helix seed
          auto hptr = art::Ptr<HelixSeed>(hseedcol_h,iseed);
          auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
          kkseedassns->addSingle(kseedptr,hptr);
          // save (unpersistable) KKTrk in the event
          ktrkcol->push_back(ktrk.release());
          //increment the counts
          if(helix_dir > 0 || (undefined_dir && index == 0)) ++nDownstream_;
          if(helix_dir < 0 || (undefined_dir && index == 1)) ++nUpstream_;
        }
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
    if(toTrackerEnds_)toTrackerEnds(ktrk,dir0);
    if(upstream_){
      TimeDir tdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
      double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
      // extrapolate through the IPA in this direction.
      bool exitsIPA = extrapolateIPA(ktrk,tdir);
      if(exitsIPA){ // if it exits out the back, extrapolate through the target
        bool exitsST = extrapolateST(ktrk,tdir);
        if(exitsST) { // if it exits out the back, extrapolate to the TSDA (DS rear absorber)
          bool hitTSDA = extrapolateTSDA(ktrk,tdir);
          // if we hit the TSDA we are done. Otherwise if we reflected, go back through the ST
          if(!hitTSDA){ // reflection upstream of the target: go back through the target
            extrapolateST(ktrk,tdir);
            if(backToTracker_){ // optionally extrapolate back through the IPA, then to the tracker entrance
              extrapolateIPA(ktrk,tdir);
              extrapolateTracker(ktrk,tdir);
            }
          }
        } else { // reflection inside the ST; extrapolate back through the IPA, then to the tracker entrance
          if(backToTracker_){
            extrapolateIPA(ktrk,tdir);
            extrapolateTracker(ktrk,tdir);
          }
        }
      } else { // reflection inside the IPA; extrapolate back through the IPA, then to the tracker entrance
        if(backToTracker_)ktrk.extrapolate(tdir,trackerFront_);
      }
      // optionally test for intersection with the OPA
      if(toOPA_)toOPA(ktrk,starttime,tdir);
    }
  }

  void LoopHelixFit::toTrackerEnds(KKTRK& ktrk,VEC3 const& dir0) const {
    TimeDir fronttdir = (dir0.Z() > 0) ? TimeDir::backwards : TimeDir::forwards;
    TimeDir backtdir = (dir0.Z() > 0) ? TimeDir::forwards : TimeDir::backwards;
    ktrk.extrapolate(fronttdir,trackerFront_);
    ktrk.extrapolate(backtdir,trackerBack_);
    // record the standard tracker intersections
    auto const& ftraj = ktrk.fitTraj();
    TimeRange frange = ftraj.range();
    double tol = trackerFront_.tolerance();
    static const SurfaceId tt_front("TT_Front");
    static const SurfaceId tt_mid("TT_Mid");
    static const SurfaceId tt_back("TT_Back");

    auto frontinter = KinKal::intersect(ftraj,*trkfrontptr_,frange,tol,fronttdir);
    if(frontinter.onsurface_)ktrk.addIntersection(tt_front,frontinter);
    auto midinter = KinKal::intersect(ftraj,*trkmidptr_,frange,tol);
    if(midinter.onsurface_)ktrk.addIntersection(tt_mid,midinter);
    auto backinter = KinKal::intersect(ftraj,*trkbackptr_,frange,tol,backtdir);
    if(backinter.onsurface_)ktrk.addIntersection(tt_back,backinter);
  }

  bool LoopHelixFit::extrapolateIPA(KKTRK& ktrk,TimeDir tdir) const {
    if(extrapIPA_.debug() > 0)std::cout << "extrapolating to IPA " << std::endl;
    // extraplate the fit through the IPA. This will add material effects for each intersection. It will continue till the
    // track exits the IPA
    extrapIPA_.reset();
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId IPASID("IPA");
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
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
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addIPAXing(ipaxingptr,tdir);
        if(extrapIPA_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapIPA_.intersection().onsurface_ && extrapIPA_.intersection().inbounds_);
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto enddir = ftraj.direction(endtime);
    if(enddir.Z() * startdir.Z() > 0.0){
      return true;
    }
    return false;
  }

  bool LoopHelixFit::extrapolateST(KKTRK& ktrk,TimeDir tdir) const {
    // extraplate the fit through the ST. This will add material effects for each foil intersection. It will continue till the
    // track exits the ST in Z
    extrapST_.reset();
    auto const& ftraj = ktrk.fitTraj();
    double starttime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto startdir = ftraj.direction(starttime);
    if(extrapST_.debug() > 0)std::cout << "extrapolating to ST " << std::endl;
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
          std::cout << " before append mom = " << reftrajptr->momentum();
        }
        ktrk.addSTXing(stxingptr,tdir);
        if(extrapST_.debug() > 0){
          auto const& newtrajptr = tdir == TimeDir::backwards ? ftraj.frontPtr() : ftraj.backPtr();
          std::cout << " after append mom = " << newtrajptr->momentum() << std::endl;
        }
      }
    } while(extrapST_.intersection().onsurface_ && extrapST_.intersection().inbounds_);
    // check if the particle exited in the same physical direction or not (reflection)
    double endtime = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto enddir = ftraj.direction(endtime);
    if(enddir.Z() * startdir.Z() > 0.0){
      return true;
    }
    return false;
  }

  bool LoopHelixFit::extrapolateTracker(KKTRK& ktrk,TimeDir tdir) const {
    if(trackerFront_.debug() > 0)std::cout << "extrapolating to Tracker " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TrackerSID("TT_Front");
    ktrk.extrapolate(tdir,trackerFront_);
    // the last piece appended should cover the necessary range
    auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
    auto trkfrontinter = KinKal::intersect(ftraj,*trkfrontptr_,ktraj.range(),trackerFront_.tolerance(),tdir);
    if(trkfrontinter.onsurface_){ // dont worry about bounds here
      ktrk.addIntersection(TrackerSID,trkfrontinter);
      return true;
    }
    return false;
  }

  bool LoopHelixFit::extrapolateTSDA(KKTRK& ktrk,TimeDir tdir) const {
    if(TSDA_.debug() > 0)std::cout << "extrapolating to TSDA " << std::endl;
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId TSDASID("TSDA");
    ktrk.extrapolate(tdir,TSDA_);
    // if we reflected we're done. Otherwize, save the TSDA intersection
    double tend = tdir == TimeDir::forwards ? ftraj.range().end() : ftraj.range().begin();
    auto epos = ftraj.position3(tend);
    bool retval = epos.Z() < TSDA_.zVal();
    if(retval){
      auto const& ktraj = tdir == TimeDir::forwards ? ftraj.back() : ftraj.front();
      auto tsdainter = KinKal::intersect(ftraj,*tsdaptr_,ktraj.range(),TSDA_.tolerance(),tdir);
      if(tsdainter.onsurface_)ktrk.addIntersection(TSDASID,tsdainter);
    }
    return retval;
  }

  void LoopHelixFit::toOPA(KKTRK& ktrk, double tstart, TimeDir tdir) const {
    auto const& ftraj = ktrk.fitTraj();
    static const SurfaceId OPASID("OPA");
    TimeRange trange = tdir == TimeDir::forwards ? TimeRange(tstart,ftraj.range().end()) : TimeRange(ftraj.range().begin(),tstart);
    auto opainter = KinKal::intersect(ftraj,*opaptr_,trange,trackerFront_.tolerance(),tdir);
    if(opainter.onsurface_ && opainter.inbounds_){ // require in bounds to say it was a physical intersection
      ktrk.addIntersection(OPASID,opainter);
    }
  }

  void LoopHelixFit::sampleFit(KKTRK const& kktrk,KalIntersectionCollection& inters) const {
    auto const& ftraj = kktrk.fitTraj();
    double tbeg = ftraj.range().begin();
    double tend = ftraj.range().end();
    static const double epsilon(1.0e-6);
    // if this helix has reflected, limit the search
    auto mom0 = ftraj.momentum3(ftraj.t0());
    if(mom0.Z() >0){ // downstream fit: skip any upstream-going segments
      for(auto const& ktraj : ftraj.pieces()){
        auto axis = ktraj->axis(ktraj->range().mid());
        if(axis.direction().Z() > 0.0 )break; // helix headed downstream
        tbeg = ktraj->range().end(); // force onto next piece
      }
    } else { // upstream fit: stop when the track reflects back downstream
      for(auto const& ktraj : ftraj.pieces()){
        auto axis = ktraj->axis(ktraj->range().mid());
        if(axis.direction().Z() > 0.0 )break; // helix no longer heading upstream
        tend = ktraj->range().begin();
      }
    }
    for(auto const& surf : sample_){
      // search for intersections with each surface within the specified time range, going forwards in time
      bool hasinter(true);
      size_t max_iter = 1000;
      size_t cur_iter = 0;
      // loop to find multiple intersections
      while(hasinter && tbeg < tend) {
        TimeRange irange(tbeg,tend);
        if (cur_iter > max_iter)
          break;
        cur_iter += 1;
        auto surfinter = KinKal::intersect(ftraj,*surf.second,irange,sampletol_);
        hasinter = surfinter.onsurface_ && ( (! sampleinbounds_) || surfinter.inbounds_ ) && ( (!sampleinrange_) || irange.inRange(surfinter.time_));
        if(hasinter) {
          // save the intersection information
          auto const& ktraj = ftraj.nearestPiece(surfinter.time_);
          inters.emplace_back(ktraj.stateEstimate(surfinter.time_),XYZVectorF(ktraj.bnom()),surf.first,surfinter);
          // update for the next intersection
          tbeg = surfinter.time_ + epsilon;// move psst existing intersection to avoid repeating
        }
      }
    }
  }

  void LoopHelixFit::endJob() {
    if(/*useHelixSlope_ &&*/ print_ > -1) printf("[LoopHelixFit::%s] Saw %i helix seeds, %i had ambiguous dz/dt slopes, accepted %i downstream and %i upstream fits\n",
                                             __func__, nSeen_, nAmbiguous_, nDownstream_, nUpstream_);
  }
}

DEFINE_ART_MODULE(mu2e::LoopHelixFit)
