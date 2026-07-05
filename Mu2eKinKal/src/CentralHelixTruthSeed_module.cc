//
// Produce a truth-seeded, extrapolated CentralHelix "track" for each reconstructed KKCH fit,
// purely from MC information, to isolate the extrapolation (B-field/material) error from the
// fit error.
//
// For every reconstructed KalSeed (CentralHelix) this module rebuilds a PIECEWISE truth helix
// with the SAME domain structure as the reco fit -- one truth helix per fit piece, each anchored
// at the true muon state (from the matched particle's MCTrajectory) at that piece's location.
// This reproduces the truth |p| -- including tracker energy loss -- at every surface, giving a
// ~zero tracker residual at all momenta; the extrapolation then carries it to the CRV. The seed
// is constrained to the truth with ParameterHits and extrapolated through the PUBLIC KinKal
// reconstruction interface. Differencing the truth-seeded CRV residual from the reco-seeded one
// gives the fit's contribution (e.g. the charge over-rotation).
//
// This is the MC-decoupled replacement for the in-module TruthSeedDiag option: it keeps all
// MC-specific code out of the CentralHelixFit reconstruction module (which nominally runs on real
// data), per D. Brown's review of PR #1869. The reco fit's piece structure and field domains are
// rebuilt from the persistent KalSeed, and the produced KalSeeds are associated with the same
// KalSeedMC as their parent reco KalSeed (cf. RegrowLoopHelix).
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
// conditions
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
// utilities
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"
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
#include <cmath>
#include <limits>

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::VEC4;
  using KinKal::MOM4;
  using KinKal::DMAT;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;

  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;

  class CentralHelixTruthSeed : public art::EDProducer {
    public:
      using KTRAJ = KinKal::CentralHelix;
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTRAJPTR = std::unique_ptr<PTRAJ>;
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
      using PARAMHIT = KinKal::ParameterHit<KTRAJ>;
      using PARAMHITPTR = std::shared_ptr<PARAMHIT>;
      using PARAMHITCOL = std::vector<PARAMHITPTR>;
      using DOMAINPTR = std::shared_ptr<KinKal::Domain>;
      using DOMAINCOL = std::set<DOMAINPTR>;
      using KKFIT = KKFit<KTRAJ>;

      struct Config {
        fhicl::Atom<int> debug{Name("debug"), Comment("Debug printout level"), 0};
        fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedCollection"), Comment("Reconstructed CentralHelix KalSeed collection to truth-seed") };
        fhicl::Atom<art::InputTag> kalSeedMCAssns {Name("KalSeedMCAssns"), Comment("Association of the reco KalSeeds to their KalSeedMC") };
        fhicl::Atom<art::InputTag> mcTrajectoryCollection {Name("MCTrajectories"), Comment("MCTrajectory map providing the true muon state along the tracker") };
        fhicl::Sequence<double> seedErrors { Name("SeedErrors"), Comment("Initial seed parameter uncertainties (rms) for d0, phi0, omega, z0, tanDip, t0") };
        fhicl::Sequence<double> truthSeedParamConstraints {
          Name("TruthSeedParameterConstraints"),
          Comment("ParameterHit RMS constraints pinning the helix to the truth: d0, phi0, omega, z0, tanDip, t0"),
          std::vector<double>{1.0e-3, 1.0e-6, 1.0e-9, 1.0e-3, 1.0e-6, 1.0e-3}
        };
        fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
        fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
        fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
        fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit CentralHelixTruthSeed(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
    private:
      int debug_;
      art::ProductToken<KalSeedCollection> kseedcol_T_;
      art::InputTag ksmca_T_;
      art::ProductToken<MCTrajectoryCollection> mctrajcol_T_;
      std::vector<double> truthSeedParamConstraints_;
      KKFIT kkfit_;
      KinKal::Config config_; // fit configuration used to build+pin the truth track
      KinKal::Config exconfig_; // extension configuration; matches the reco track's config at extrapolation
      DMAT seedcov_; // seed covariance matrix
      std::unique_ptr<KKBField> kkbf_;
      std::unique_ptr<KKExtrap> extrap_;
  };

  CentralHelixTruthSeed::CentralHelixTruthSeed(const Parameters& settings) : art::EDProducer(settings),
    debug_(settings().debug()),
    kseedcol_T_(consumes<KalSeedCollection>(settings().kalSeedCollection())),
    ksmca_T_(settings().kalSeedMCAssns()),
    mctrajcol_T_(consumes<MCTrajectoryCollection>(settings().mcTrajectoryCollection())),
    truthSeedParamConstraints_(settings().truthSeedParamConstraints()),
    kkfit_(settings().kkfitSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings()))
  {
    consumes<KalSeedMCAssns>(ksmca_T_);
    produces<KKTRKCOL>();
    produces<KalSeedCollection>();
    produces<KalSeedMCAssns>();
    if(settings().extrapSettings()) extrap_ = std::make_unique<KKExtrap>(*settings().extrapSettings());
    if(truthSeedParamConstraints_.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: TruthSeedParameterConstraints must have "
        << KinKal::NParams() << " entries" << std::endl;
    for(auto const& sigma : truthSeedParamConstraints_)
      if(sigma <= 0.0) throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: TruthSeedParameterConstraints entries must be positive" << std::endl;
    auto const& seederrors = settings().seedErrors();
    if(seederrors.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: SeedErrors must have " << KinKal::NParams() << " entries" << std::endl;
    for(size_t ipar = 0; ipar < seederrors.size(); ++ipar) seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
  }

  void CentralHelixTruthSeed::beginRun(art::Run& run) {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::make_unique<KKBField>(*bfmgr,*det);
  }

  void CentralHelixTruthSeed::produce(art::Event& event) {
    GeomHandle<Calorimeter> calo_h;
    GeomHandle<mu2e::Tracker> nominalTracker_h;
    GeomHandle<DetectorSystem> det;
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();

    auto kseed_H = event.getValidHandle<KalSeedCollection>(kseedcol_T_);
    auto const& kseedcol = *kseed_H;
    auto ksmca_H = event.getValidHandle<KalSeedMCAssns>(ksmca_T_);
    auto const& ksmca = *ksmca_H;
    auto const* mctrajs = event.getValidHandle<MCTrajectoryCollection>(mctrajcol_T_).product();
    if(ksmca.size() != kseedcol.size())
      throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: KalSeedMCAssns size " << ksmca.size()
        << " does not match KalSeedCollection size " << kseedcol.size() << std::endl;

    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    auto ktrkcol = std::make_unique<KKTRKCOL>();
    auto tsseedcol = std::make_unique<KalSeedCollection>();
    auto tsmca = std::make_unique<KalSeedMCAssns>();

    // true muon state (pos, mom, time) from the MC trajectory point nearest a given detector-frame location.
    // NB: do NOT dereference the trajectory's SimParticle Ptr -- cosmic resampling/mixing leaves it dangling
    // (ProductNotFound); the trajectory POINTS (pos/time/KE) are valid.
    auto trueStateAt = [&](VEC3 const& at, double mass, VEC3& tpos, VEC3& tmom, double& ttime)->bool{
      double bestd2 = std::numeric_limits<double>::max(); bool found = false;
      for(auto const& tjpair : *mctrajs){
        auto const& pts = tjpair.second.points();
        if(pts.size() < 2) continue;
        for(size_t i = 0; i < pts.size(); ++i){
          auto pdh = det->toDetector(pts[i].pos());
          VEC3 pd(pdh.x(),pdh.y(),pdh.z());
          double d2 = (pd - at).Mag2();
          if(d2 >= bestd2) continue;
          size_t j = (i + 1 < pts.size()) ? i + 1 : i - 1; // local tangent for the direction
          auto qdh = det->toDetector(pts[j].pos());
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
    };

    for(size_t iseed = 0; iseed < kseedcol.size(); ++iseed) {
      auto const& kseed = kseedcol[iseed];
      if(!kseed.centralHelixFit())
        throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: passed KalSeed from non-CentralHelix fit" << std::endl;
      auto const& ksmcai = ksmca[iseed];
      auto origksp = art::Ptr<KalSeed>(kseed_H,iseed);
      if(ksmcai.first != origksp)
        throw cet::exception("RECO") << "mu2e::CentralHelixTruthSeed: KalSeedMCAssns is not parallel to the KalSeedCollection" << std::endl;
      auto const& mcseedp = ksmcai.second;
      double mass = ptable->particle(kseed.particle()).mass();

      // reconstruct the reco fit trajectory (piecewise) from the persistent seed
      auto recotraj = kseed.centralHelixFitTrajectory();
      int tcharge = (recotraj->nearestPiece(recotraj->t0()).charge() > 0.0) ? 1 : -1;

      // Rebuild the truth seed with the SAME domain structure as the fit: one truth helix per fit piece,
      // anchored at the true muon state at that piece's location. This reproduces the truth |p| (tracker
      // energy loss included) at every surface.
      PTRAJ truthpt;
      std::vector<double> truthHitTimes;
      for(auto const& fpieceptr : recotraj->pieces()){
        auto const& fpiece = *fpieceptr;
        auto const htime = fpiece.range().mid();
        VEC3 tpos, tmom; double ttime = 0.0;
        if(!trueStateAt(fpiece.position3(htime), mass, tpos, tmom, ttime)) continue;
        auto tbnom = kkbf_->fieldVect(tpos);
        KTRAJ tpiece(VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime),
                     MOM4(tmom.X(),tmom.Y(),tmom.Z(),mass),
                     tcharge, tbnom.Z(), fpiece.range());
        tpiece.params() = KinKal::Parameters(tpiece.params().parameters(), seedcov_);
        truthpt.append(tpiece);
        truthHitTimes.emplace_back(htime);
      }
      if(truthHitTimes.empty()) continue;

      // constrain each truth piece with a ParameterHit at its reference time (public reco interface).
      PARAMHITCOL truthParamHits;
      PARAMHIT::PMASK truthMask; truthMask.fill(true);
      auto addTruthHit = [&](double hitTime, double covarianceScale){
        auto cparams = truthpt.nearestPiece(hitTime).params();
        for(size_t ipar = 0; ipar < KinKal::NParams(); ++ipar){
          for(size_t jpar = 0; jpar < KinKal::NParams(); ++jpar) cparams.covariance()[ipar][jpar] = 0.0;
          auto const sigma = truthSeedParamConstraints_.at(ipar);
          cparams.covariance()[ipar][ipar] = covarianceScale*sigma*sigma;
        }
        truthParamHits.push_back(std::make_shared<PARAMHIT>(hitTime, truthpt, cparams, truthMask));
      };
      // a single piece gives a zero-length fit range (lowNDOF); split it into two hits a sliver apart,
      // with looser (x2) constraint.
      if(truthHitTimes.size() == 1){
        constexpr double truthHitDt = 1.0e-3;
        addTruthHit(truthHitTimes.front() - 0.5*truthHitDt, 2.0);
        addTruthHit(truthHitTimes.front() + 0.5*truthHitDt, 2.0);
      } else {
        for(auto const hitTime : truthHitTimes) addTruthHit(hitTime, 1.0);
      }
      // rebuild the fit's field domains from the seed so the truth fit/extrapolation see the same BField
      // model (same recipe as KKFit::regrowComponents)
      DOMAINCOL truthDomains;
      auto const& dbounds = kseed.domainBounds();
      for(size_t idb = 0; idb + 1 < dbounds.size(); ++idb){
        double tstart = dbounds[idb];
        double trange = dbounds[idb+1] - tstart;
        double tmid = tstart + 0.5*trange;
        truthDomains.emplace(std::make_shared<KinKal::Domain>(tstart,trange,recotraj->nearestPiece(tmid).bnom()));
      }

      PTRAJPTR truthTraj = std::make_unique<PTRAJ>(truthpt);
      KKSTRAWHITCOL nostrawhits; KKSTRAWXINGCOL nostrawxings; KKCALOHITCOL nocalohits;
      auto truthConfig = config_;
      truthConfig.maxniter_ = 1;
      truthConfig.divdchisq_ = std::numeric_limits<double>::max();
      truthConfig.pdchisq_ = std::numeric_limits<double>::max();
      truthConfig.divgap_ = std::numeric_limits<double>::max();
      // Extrapolation must be IDENTICAL to the reco track's. The reco fit is extended with the
      // ExtensionSettings config before extrapolating, so its extrapolation field handling is governed by
      // exconfig_.bfcorr_ (resample BField vs freeze BNom) and exconfig_.tol_ (BField-correction tolerance).
      // Carry those two fields onto the truth track's config so the extrapolation machinery is identical --
      // WITHOUT re-running the full extension schedule, which diverges a hit-less truth-pinned fit. The
      // pinning itself keeps the (proven) seed-fit schedule; only the extrapolation-relevant fields change.
      truthConfig.bfcorr_ = exconfig_.bfcorr_;
      truthConfig.tol_    = exconfig_.tol_;
      auto ktrk_truth = std::make_unique<KKTRK>(truthConfig,*kkbf_,kseed.particle(),truthTraj,
          nostrawhits,nostrawxings,nocalohits,truthParamHits,truthDomains);
      if(ktrk_truth->fitStatus().usable()){
        if(extrap_) extrap_->extrapolate(*ktrk_truth);
        // Sample the (extrapolated) truth track at the configured surfaces (tracker ends + CRV sectors),
        // exactly as the reco does (via RegrowKalSeed's kkfit_.sampleFit). The truth track has no straw hits,
        // so without this the KalSeed's detector time range is undefined ("tracker intersections missing")
        // and the CRV crossings are not recorded -- leaving trktruthsegs empty at the CRV.
        kkfit_.sampleFit(*ktrk_truth);
        TrkFitFlag tflag(kseed.status()); tflag.merge(TrkFitFlag::KKCentralHelix); tflag.merge(TrkFitFlag::FitOK);
        auto tsseed = kkfit_.createSeed(*ktrk_truth,tflag,*calo_h,*nominalTracker_h);
        // carry the reco fit's hit metadata onto the truth KalSeed (bookkeeping for EventNtuple consumers)
        tsseed._hits = kseed.hits();
        tsseed._hitcalibs = kseed.hitCalibInfos();
        tsseed._chit = kseed.caloHit();
        tsseedcol->push_back(std::move(tsseed));
        auto tsseedp = art::Ptr<KalSeed>(KalSeedCollectionPID,tsseedcol->size()-1,KalSeedCollectionGetter);
        tsmca->addSingle(tsseedp,mcseedp);
        ktrkcol->push_back(ktrk_truth.release());
      } else if(debug_ > 0){
        std::cout << "CentralHelixTruthSeed ParameterHit fit unusable: " << ktrk_truth->fitStatus() << std::endl;
      }
    }
    if(debug_ > 0) std::cout << "CentralHelixTruthSeed made " << tsseedcol->size()
      << " truth-seeded tracks from " << kseedcol.size() << " reco KalSeeds" << std::endl;
    event.put(std::move(ktrkcol));
    event.put(std::move(tsseedcol));
    event.put(std::move(tsmca));
  }
}
DEFINE_ART_MODULE(mu2e::CentralHelixTruthSeed)
