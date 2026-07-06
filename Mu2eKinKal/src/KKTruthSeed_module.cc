//
// Produce a truth-seeded, extrapolated "track" for each reconstructed KalSeed, purely from MC
// information, to isolate the extrapolation (B-field / material / scattering) error from the fit error.
//
// A single module handles every supported fit type (LoopHelix, CentralHelix, KinematicLine): for each
// reconstructed KalSeed it dispatches on the seed's fit type to a templated worker that rebuilds a
// parallel track seeded with the true muon state -- sampled from the matched
// particle's MCTrajectory at each reco-fit-piece location -- constrains it to that truth with
// ParameterHits, and extrapolates it through the public KinKal reconstruction interface using the
// identical fit/extension/extrapolation configuration as the reco track. Differencing the truth-seeded
// CRV residual from the reco-seeded one gives the fit's contribution to the residual.
//
// Original author: R. Mina (UVA) 2026, following D. Brown's RegrowKalSeed pattern
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
// conditions / geometry
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKExtrap.hh"
// C++
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <limits>
#include <cmath>
#include <type_traits>

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::VEC4;
  using KinKal::MOM4;
  using KinKal::DMAT;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;

  class KKTruthSeed : public art::EDProducer {
    public:
      struct Config {
        fhicl::Atom<int> debug{Name("debug"), Comment("Debug printout level"), 0};
        fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedCollection"), Comment("Reconstructed KalSeed collection to truth-seed") };
        fhicl::Atom<art::InputTag> kalSeedMCAssns {Name("KalSeedMCAssns"), Comment("Association of the reco KalSeeds to their KalSeedMC (the truth match)") };
        fhicl::Atom<art::InputTag> mcTrajectoryCollection {Name("MCTrajectories"), Comment("MCTrajectory map providing the true muon state along the track") };
        fhicl::Sequence<double> seedErrors { Name("SeedErrors"), Comment("Initial seed parameter uncertainties (rms)") };
        fhicl::Sequence<double> truthSeedParamConstraints {
          Name("TruthSeedParameterConstraints"),
          Comment("ParameterHit RMS constraints pinning the seed to the truth (one per fit parameter)"),
          std::vector<double>{1.0e-3, 1.0e-6, 1.0e-9, 1.0e-3, 1.0e-6, 1.0e-3}
        };
        fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
        fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
        fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
        fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit KKTruthSeed(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;

    private:
      // true muon state (pos, mom, time) from the MC trajectory point nearest a detector-frame location.
      // NB: do NOT dereference the trajectory's SimParticle Ptr -- cosmic resampling/mixing leaves it
      // dangling (ProductNotFound); the trajectory POINTS (pos/time/KE) are valid.
      bool trueStateAt(MCTrajectoryCollection const& mctrajs, DetectorSystem const& det,
                       VEC3 const& at, double mass, VEC3& tpos, VEC3& tmom, double& ttime) const;

      // templated worker: build the truth-seeded track for one reco KalSeed of type KTRAJ (extrapolating
      // it if configured), returning its truth KalSeed (or nullopt if it can't be built / the fit is unusable).
      template <class KTRAJ>
      std::optional<KalSeed> createTruthSeed(KinKal::ParticleTrajectory<KTRAJ> const& recotraj,
          KalSeed const& recoseed, KKFit<KTRAJ>& kkfit, MCTrajectoryCollection const& mctrajs,
          DetectorSystem const& det, ParticleDataList const& ptable,
          Calorimeter const& calo, Tracker const& nominalTracker) const;

      int debug_;
      art::ProductToken<KalSeedCollection> kseedcol_T_;
      art::InputTag ksmca_T_;
      art::ProductToken<MCTrajectoryCollection> mctraj_T_;
      std::vector<double> truthSeedParamConstraints_;
      KinKal::Config config_;   // fit configuration used to build+pin the truth track
      KinKal::Config exconfig_; // extension configuration; supplies the extrapolation-relevant fields
      DMAT seedcov_;            // seed covariance matrix
      std::unique_ptr<KKBField> kkbf_;
      std::unique_ptr<KKExtrap> extrap_;
      KKFit<KinKal::CentralHelix>  chfit_;
      KKFit<KinKal::KinematicLine> klfit_;
      KKFit<KinKal::LoopHelix>     lhfit_;
  };

  KKTruthSeed::KKTruthSeed(const Parameters& settings) : art::EDProducer(settings),
    debug_(settings().debug()),
    kseedcol_T_(consumes<KalSeedCollection>(settings().kalSeedCollection())),
    ksmca_T_(settings().kalSeedMCAssns()),
    mctraj_T_(consumes<MCTrajectoryCollection>(settings().mcTrajectoryCollection())),
    truthSeedParamConstraints_(settings().truthSeedParamConstraints()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings())),
    chfit_(settings().kkfitSettings()),
    klfit_(settings().kkfitSettings()),
    lhfit_(settings().kkfitSettings())
  {
    consumes<KalSeedMCAssns>(ksmca_T_);
    produces<KalSeedCollection>();
    produces<KalSeedMCAssns>();
    if(settings().extrapSettings()) extrap_ = std::make_unique<KKExtrap>(*settings().extrapSettings());
    if(truthSeedParamConstraints_.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::KKTruthSeed: TruthSeedParameterConstraints must have "
        << KinKal::NParams() << " entries" << std::endl;
    for(auto const& sigma : truthSeedParamConstraints_)
      if(sigma <= 0.0) throw cet::exception("RECO") << "mu2e::KKTruthSeed: TruthSeedParameterConstraints entries must be positive" << std::endl;
    auto const& seederrors = settings().seedErrors();
    if(seederrors.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::KKTruthSeed: SeedErrors must have " << KinKal::NParams() << " entries" << std::endl;
    for(size_t ipar = 0; ipar < seederrors.size(); ++ipar) seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
  }

  void KKTruthSeed::beginRun(art::Run& run) {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::make_unique<KKBField>(*bfmgr,*det);
  }

  bool KKTruthSeed::trueStateAt(MCTrajectoryCollection const& mctrajs, DetectorSystem const& det,
      VEC3 const& at, double mass, VEC3& tpos, VEC3& tmom, double& ttime) const {
    double bestd2 = std::numeric_limits<double>::max(); bool found = false;
    for(auto const& tjpair : mctrajs){
      auto const& pts = tjpair.second.points();
      if(pts.size() < 2) continue;
      for(size_t i = 0; i < pts.size(); ++i){
        auto pdh = det.toDetector(pts[i].pos());
        VEC3 pd(pdh.x(),pdh.y(),pdh.z());
        double d2 = (pd - at).Mag2();
        if(d2 >= bestd2) continue;
        size_t j = (i + 1 < pts.size()) ? i + 1 : i - 1; // local tangent for the direction
        auto qdh = det.toDetector(pts[j].pos());
        VEC3 qd(qdh.x(),qdh.y(),qdh.z());
        VEC3 dir = (j > i) ? (qd - pd) : (pd - qd);
        double dn = std::sqrt(dir.Mag2());
        if(dn <= 0.0) continue;
        dir /= dn;
        double ke = pts[i].kineticEnergy();
        double pmag = std::sqrt(ke*ke + 2.0*ke*mass); // |p| from kinetic energy
        bestd2 = d2; found = true; tpos = pd; tmom = pmag*dir; ttime = pts[i].t();
      }
    }
    return found;
  }

  template <class KTRAJ>
  std::optional<KalSeed> KKTruthSeed::createTruthSeed(KinKal::ParticleTrajectory<KTRAJ> const& recotraj,
      KalSeed const& recoseed, KKFit<KTRAJ>& kkfit, MCTrajectoryCollection const& mctrajs,
      DetectorSystem const& det, ParticleDataList const& ptable,
      Calorimeter const& calo, Tracker const& nominalTracker) const {
    using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
    using KKTRK = KKTrack<KTRAJ>;
    using KKSTRAWHITCOL = std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>;
    using KKSTRAWXINGCOL = std::vector<std::shared_ptr<KKStrawXing<KTRAJ>>>;
    using KKCALOHITCOL = std::vector<std::shared_ptr<KKCaloHit<KTRAJ>>>;
    using PARAMHIT = KinKal::ParameterHit<KTRAJ>;
    using PARAMHITCOL = std::vector<std::shared_ptr<PARAMHIT>>;
    using DOMAINCOL = std::set<std::shared_ptr<KinKal::Domain>>;
    // KinematicLine is the only field-FREE fit: it is a single trajectory piece with no BField domains.
    constexpr bool isLine = std::is_same_v<KTRAJ,KinKal::KinematicLine>;

    double mass = ptable.particle(recoseed.particle()).mass();
    int tcharge = (recotraj.nearestPiece(recotraj.t0()).charge() > 0.0) ? 1 : -1;

    // Rebuild the truth seed with the SAME piece structure as the fit: one truth piece per fit piece,
    // anchored at the true muon state at that piece's location. For a helical (field-on) fit this
    // reproduces the truth |p| (tracker energy loss included) at every surface; for a single-piece line
    // (field-off) it is just the true tracker state.
    PTRAJ truthpt;
    std::vector<double> truthHitTimes;
    for(auto const& fpieceptr : recotraj.pieces()){
      auto const& fpiece = *fpieceptr;
      auto const htime = fpiece.range().mid();
      VEC3 tpos, tmom; double ttime = 0.0;
      if(!trueStateAt(mctrajs, det, fpiece.position3(htime), mass, tpos, tmom, ttime)) continue;
      // the KTRAJ nominal BField: the local field for a helix; a fixed non-zero z for a field-free line
      // (KinematicLine requires a non-zero bnom for interface consistency).
      KTRAJ tpiece = [&]{
        if constexpr (isLine){
          return KTRAJ(VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime), MOM4(tmom.X(),tmom.Y(),tmom.Z(),mass),
                       tcharge, VEC3(0.0,0.0,0.001), fpiece.range());
        } else {
          return KTRAJ(VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime), MOM4(tmom.X(),tmom.Y(),tmom.Z(),mass),
                       tcharge, kkbf_->fieldVect(tpos).Z(), fpiece.range());
        }
      }();
      tpiece.params() = KinKal::Parameters(tpiece.params().parameters(), seedcov_);
      truthpt.append(tpiece);
      truthHitTimes.emplace_back(htime);
    }
    if(truthHitTimes.empty()) return std::nullopt;

    // constrain each truth piece with a ParameterHit at its (in-range) reference time -- the public reco
    // interface (no fitTraj substitution). Piece midpoints are always in-range, so no time clamping.
    PARAMHITCOL truthParamHits;
    typename PARAMHIT::PMASK truthMask; truthMask.fill(true);
    auto addTruthHit = [&](double hitTime, double covarianceScale){
      auto cparams = truthpt.nearestPiece(hitTime).params();
      for(size_t ipar = 0; ipar < KinKal::NParams(); ++ipar){
        for(size_t jpar = 0; jpar < KinKal::NParams(); ++jpar) cparams.covariance()[ipar][jpar] = 0.0;
        auto const sigma = truthSeedParamConstraints_.at(ipar);
        cparams.covariance()[ipar][ipar] = covarianceScale*sigma*sigma;
      }
      truthParamHits.push_back(std::make_shared<PARAMHIT>(hitTime, truthpt, cparams, truthMask));
    };
    if(truthHitTimes.size() == 1){
      // Rare degenerate case: recotraj is normally MANY pieces (one per material/BField crossing), so the
      // per-piece branch below runs. But a lone piece gives a single ParameterHit == a zero-length (lowNDOF)
      // fit range, so pin two hits a sliver apart instead.
      constexpr double truthHitDt = 1.0e-3;
      addTruthHit(truthHitTimes.front() - 0.5*truthHitDt, 2.0);
      addTruthHit(truthHitTimes.front() + 0.5*truthHitDt, 2.0);
    } else {
      // pin each trajectory piece at its reference-time midpoint.
      for(auto const hitTime : truthHitTimes) addTruthHit(hitTime, 1.0);
    }

    // rebuild the fit's field domains from the seed so the truth fit/extrapolation see the same BField
    // model. Field-on fits only: a field-free line (KinematicLine) has no BField domains.
    DOMAINCOL truthDomains;
    if constexpr(!isLine){
      auto const& dbounds = recoseed.domainBounds();
      for(size_t idb = 0; idb + 1 < dbounds.size(); ++idb){
        double tstart = dbounds[idb];
        double trange = dbounds[idb+1] - tstart;
        double tmid = tstart + 0.5*trange;
        truthDomains.emplace(std::make_shared<KinKal::Domain>(tstart,trange,recotraj.nearestPiece(tmid).bnom()));
      }
    }

    std::unique_ptr<PTRAJ> truthTraj = std::make_unique<PTRAJ>(truthpt);
    KKSTRAWHITCOL nostrawhits; KKSTRAWXINGCOL nostrawxings; KKCALOHITCOL nocalohits;
    // one iteration: the ParameterHits already pin the params.
    auto truthConfig = config_;
    truthConfig.maxniter_ = 1;
    truthConfig.divdchisq_ = std::numeric_limits<double>::max();
    truthConfig.pdchisq_ = std::numeric_limits<double>::max();
    truthConfig.divgap_ = std::numeric_limits<double>::max();
    // Extrapolation should be identical to the reco track's. The reco fit is extended with the
    // ExtensionSettings config before extrapolating, so its extrapolation field handling is governed by
    // exconfig_.bfcorr_ (resample BField vs freeze BNom) and exconfig_.tol_ (BField-correction tolerance).
    // Carry those two fields onto the truth track's config so the extrapolation machinery is identical --
    // without re-running the full extension schedule (which diverges a hit-less truth-pinned fit).
    truthConfig.bfcorr_ = exconfig_.bfcorr_;
    truthConfig.tol_    = exconfig_.tol_;
    auto ktrk_truth = std::make_unique<KKTRK>(truthConfig,*kkbf_,recoseed.particle(),truthTraj,
        nostrawhits,nostrawxings,nocalohits,truthParamHits,truthDomains);
    if(!ktrk_truth->fitStatus().usable()){
      if(debug_ > 0) std::cout << "KKTruthSeed ParameterHit fit unusable: " << ktrk_truth->fitStatus() << std::endl;
      return std::nullopt;
    }
    if(extrap_) extrap_->extrapolate(*ktrk_truth);
    // Sample the (extrapolated) truth track at the configured surfaces (tracker ends + CRV sectors), as
    // the reco does. The truth track has no straw hits, so without this the KalSeed's detector time range
    // is undefined and the CRV crossings are not recorded -- leaving trktruthsegs empty at the CRV.
    kkfit.sampleFit(*ktrk_truth);
    auto tsseed = kkfit.createSeed(*ktrk_truth,recoseed.status(),calo,nominalTracker);
    // carry the reco fit's hit metadata onto the truth KalSeed (the truth-constrained fit has no detector
    // hits, but EventNtuple consumers use the hit list for bookkeeping)
    tsseed._hits = recoseed.hits();
    tsseed._hitcalibs = recoseed.hitCalibInfos();
    tsseed._chit = recoseed.caloHit();
    return tsseed;
  }

  void KKTruthSeed::produce(art::Event& event) {
    GeomHandle<Calorimeter> calo_h;
    GeomHandle<mu2e::Tracker> nominalTracker_h;
    GeomHandle<DetectorSystem> det;
    auto const& ptable = *GlobalConstantsHandle<ParticleDataList>();

    auto kseed_H = event.getValidHandle<KalSeedCollection>(kseedcol_T_);
    auto const& kseedcol = *kseed_H;
    auto ksmca_H = event.getValidHandle<KalSeedMCAssns>(ksmca_T_);
    auto const& ksmca = *ksmca_H;
    auto const& mctrajs = *event.getValidHandle<MCTrajectoryCollection>(mctraj_T_);
    if(ksmca.size() != kseedcol.size())
      throw cet::exception("RECO") << "mu2e::KKTruthSeed: KalSeedMCAssns size " << ksmca.size()
        << " does not match KalSeedCollection size " << kseedcol.size() << std::endl;

    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    auto tsseedcol = std::make_unique<KalSeedCollection>();
    auto tsmca = std::make_unique<KalSeedMCAssns>();

    for(size_t iseed = 0; iseed < kseedcol.size(); ++iseed) {
      auto const& kseed = kseedcol[iseed];
      // the KalSeedMC associated to this reco KalSeed (the truth match, and the association we propagate)
      auto const& ksmcai = ksmca[iseed];
      auto origksp = art::Ptr<KalSeed>(kseed_H,iseed);
      if(ksmcai.first != origksp)
        throw cet::exception("RECO") << "mu2e::KKTruthSeed: KalSeedMCAssns is not parallel to the KalSeedCollection" << std::endl;
      auto const& mcseedp = ksmcai.second;

      // dispatch on the reco fit type (the runtime type check): build the matching KTRAJ truth seed
      std::optional<KalSeed> tsseed;
      if(kseed.loopHelixFit()){
        auto recotraj = kseed.loopHelixFitTrajectory();
        tsseed = createTruthSeed<KinKal::LoopHelix>(*recotraj,kseed,lhfit_,mctrajs,*det,ptable,*calo_h,*nominalTracker_h);
      } else if(kseed.centralHelixFit()){
        auto recotraj = kseed.centralHelixFitTrajectory();
        tsseed = createTruthSeed<KinKal::CentralHelix>(*recotraj,kseed,chfit_,mctrajs,*det,ptable,*calo_h,*nominalTracker_h);
      } else if(kseed.kinematicLineFit()){
        auto recotraj = kseed.kinematicLineFitTrajectory();
        tsseed = createTruthSeed<KinKal::KinematicLine>(*recotraj,kseed,klfit_,mctrajs,*det,ptable,*calo_h,*nominalTracker_h);
      } else {
        throw cet::exception("RECO") << "mu2e::KKTruthSeed: KalSeed is not a LoopHelix, CentralHelix, or KinematicLine fit" << std::endl;
      }
      if(!tsseed) continue;

      tsseedcol->push_back(std::move(*tsseed));
      // associate the truth-seeded KalSeed to the SAME KalSeedMC as its reco parent
      auto tsseedp = art::Ptr<KalSeed>(KalSeedCollectionPID,tsseedcol->size()-1,KalSeedCollectionGetter);
      tsmca->addSingle(tsseedp,mcseedp);
    }
    if(debug_ > 0) std::cout << "KKTruthSeed made " << tsseedcol->size()
      << " truth-seeded tracks from " << kseedcol.size() << " reco KalSeeds" << std::endl;
    event.put(std::move(tsseedcol));
    event.put(std::move(tsmca));
  }
}
DEFINE_ART_MODULE(mu2e::KKTruthSeed)
