//
// Produce a truth-seeded, extrapolated KinematicLine "track" for each reconstructed
// KKLine fit, purely from MC information, to isolate the extrapolation (material/scattering)
// error from the fit error.
//
// For every reconstructed KalSeed (KKLine) this module builds a parallel KinematicLine
// seeded with the TRUE muon state at the tracker -- taken from the matched KalSeedMC's
// VD steps -- constrains it to that truth with ParameterHits, and extrapolates it through
// the PUBLIC KinKal reconstruction interface. Differencing the truth-seeded CRV residual
// from the reco-seeded one gives the fit's contribution to the residual.
//
// This is the MC-decoupled replacement for the in-module TruthSeedDiag option: it keeps
// all MC-specific code out of the KinematicLineFit reconstruction module (which nominally
// runs on real data), per D. Brown's review of PR #1869. The produced KalSeeds are
// associated with the same KalSeedMC as their parent reco KalSeed, so downstream consumers
// can relate truth-seeded and reco-seeded tracks via the MC match (cf. RegrowKinematicLine).
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
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
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
#include <limits>

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::DMAT;
  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;

  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;

  class KinematicLineTruthSeed : public art::EDProducer {
    public:
      using KTRAJ = KinKal::KinematicLine;
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
        fhicl::Atom<art::InputTag> kalSeedCollection {Name("KalSeedCollection"), Comment("Reconstructed KKLine KalSeed collection to truth-seed") };
        fhicl::Atom<art::InputTag> kalSeedMCAssns {Name("KalSeedMCAssns"), Comment("Association of the reco KalSeeds to their KalSeedMC (source of the true muon VD state)") };
        fhicl::Sequence<double> seedErrors { Name("SeedErrors"), Comment("Initial seed parameter uncertainties (rms) for d0, phi0, z0, cost, t0, mom") };
        fhicl::Sequence<double> truthSeedParamConstraints {
          Name("TruthSeedParameterConstraints"),
          Comment("ParameterHit RMS constraints pinning the line to the truth: d0, phi0, z0, cost, t0, mom"),
          std::vector<double>{1.0e-3, 1.0e-6, 1.0e-3, 1.0e-6, 1.0e-3, 1.0e-3}
        };
        fhicl::Table<KKFitConfig> kkfitSettings { Name("KKFitSettings") };
        fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
        fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
        fhicl::OptionalTable<KKExtrapConfig> extrapSettings { Name("ExtrapolationSettings") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit KinematicLineTruthSeed(const Parameters& settings);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
    private:
      int debug_;
      art::ProductToken<KalSeedCollection> kseedcol_T_;
      art::InputTag ksmca_T_;
      std::vector<double> truthSeedParamConstraints_;
      KKFIT kkfit_;
      KinKal::Config config_; // fit configuration used to build+pin the truth track
      KinKal::Config exconfig_; // extension configuration; matches the reco track's config at extrapolation
      DMAT seedcov_; // seed covariance matrix
      std::unique_ptr<KKBField> kkbf_;
      std::unique_ptr<KKExtrap> extrap_;
      double mass_ = 0.0;
  };

  KinematicLineTruthSeed::KinematicLineTruthSeed(const Parameters& settings) : art::EDProducer(settings),
    debug_(settings().debug()),
    kseedcol_T_(consumes<KalSeedCollection>(settings().kalSeedCollection())),
    ksmca_T_(settings().kalSeedMCAssns()),
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
    // validate the truth constraints
    if(truthSeedParamConstraints_.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: TruthSeedParameterConstraints must have "
        << KinKal::NParams() << " entries" << std::endl;
    for(auto const& sigma : truthSeedParamConstraints_)
      if(sigma <= 0.0) throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: TruthSeedParameterConstraints entries must be positive" << std::endl;
    // build the initial seed covariance
    auto const& seederrors = settings().seedErrors();
    if(seederrors.size() != KinKal::NParams())
      throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: SeedErrors must have " << KinKal::NParams() << " entries" << std::endl;
    for(size_t ipar = 0; ipar < seederrors.size(); ++ipar) seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
  }

  void KinematicLineTruthSeed::beginRun(art::Run& run) {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::make_unique<KKBField>(*bfmgr,*det);
  }

  void KinematicLineTruthSeed::produce(art::Event& event) {
    GeomHandle<Calorimeter> calo_h;
    GeomHandle<mu2e::Tracker> nominalTracker_h;
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();

    // inputs: reco KalSeeds and their MC match
    auto kseed_H = event.getValidHandle<KalSeedCollection>(kseedcol_T_);
    auto const& kseedcol = *kseed_H;
    auto ksmca_H = event.getValidHandle<KalSeedMCAssns>(ksmca_T_);
    auto const& ksmca = *ksmca_H;
    if(ksmca.size() != kseedcol.size())
      throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: KalSeedMCAssns size " << ksmca.size()
        << " does not match KalSeedCollection size " << kseedcol.size() << std::endl;

    // outputs
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    auto ktrkcol = std::make_unique<KKTRKCOL>();
    auto tsseedcol = std::make_unique<KalSeedCollection>();
    auto tsmca = std::make_unique<KalSeedMCAssns>();

    for(size_t iseed = 0; iseed < kseedcol.size(); ++iseed) {
      auto const& kseed = kseedcol[iseed];
      if(!kseed.kinematicLineFit())
        throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: passed KalSeed from non-KinematicLine fit" << std::endl;
      // the KalSeedMC associated to this reco KalSeed (source of the true muon state)
      auto const& ksmcai = ksmca[iseed];
      auto origksp = art::Ptr<KalSeed>(kseed_H,iseed);
      if(ksmcai.first != origksp)
        throw cet::exception("RECO") << "mu2e::KinematicLineTruthSeed: KalSeedMCAssns is not parallel to the KalSeedCollection" << std::endl;
      auto const& mcseedp = ksmcai.second;
      auto const& mcseed = *mcseedp;
      if(mcseed.simParticles().empty() || mcseed.vdSteps().empty()) continue;

      // the reco fit trajectory gives the reference point (t0) and time range to seed against
      auto recotraj = kseed.kinematicLineFitTrajectory();
      double tt0 = recotraj->t0();
      auto refpos = recotraj->position3(tt0);

      // pick the matched particle's VD crossing nearest the fit t0 point as the true tracker state
      int tpdg = mcseed.simParticle(0)._pdg;
      double tmass = ptable->particle(static_cast<PDGCode::type>(tpdg)).mass();
      int tcharge = (tpdg > 0) ? -1 : 1; // mu-(+13) -> -1, mu+(-13) -> +1
      double bestd2 = std::numeric_limits<double>::max();
      VEC3 tpos, tmom; double ttime = 0.0; bool found = false;
      for(auto const& vds : mcseed.vdSteps()){ // vdStep positions are already in DETECTOR coordinates
        VEC3 dpos(vds._pos.x(),vds._pos.y(),vds._pos.z());
        double d2 = (dpos-refpos).Mag2();
        if(d2 < bestd2){
          bestd2 = d2; found = true; ttime = vds._time;
          tpos = dpos; tmom = VEC3(vds._mom.x(),vds._mom.y(),vds._mom.z());
        }
      }
      if(!found) continue;

      try {
        // field-off straight line: same non-zero nominal bfield (VEC3) as the KKLine seed traj
        VEC3 tbnom(0.0,0.0,0.001);
        KTRAJ truthtraj(KinKal::VEC4(tpos.X(),tpos.Y(),tpos.Z(),ttime),
                        KinKal::MOM4(tmom.X(),tmom.Y(),tmom.Z(),tmass),
                        tcharge, tbnom, recotraj->range());
        truthtraj.params() = KinKal::Parameters(truthtraj.params().parameters(), seedcov_);
        PTRAJ truthpt(truthtraj);
        // constrain the line to the truth params via ParameterHits (public reco interface). A single
        // piece is a zero-length fit range (lowNDOF), so pin TWO hits a sliver apart.
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
        constexpr double truthHitDt = 1.0e-3;
        // Pin the constraints at an IN-RANGE time. The KKLine fit is field-off, so there are no field
        // domains and Track::fit's detrange collapses to just these hit times. The VD crossing time
        // `ttime` can fall outside the fit trajectory's time range; pinning the hits there would leave
        // the single truth piece not overlapping detrange, so PiecewiseTrajectory::setRange(trim) would
        // erase every piece and dereference the empty deque (SIGSEGV). The single-piece KinematicLine
        // params are time-independent, so clamping the constraint time into the fit range pins exactly
        // the same truth params.
        auto const& frange = recotraj->range();
        double tctr = std::max(frange.begin() + truthHitDt, std::min(ttime, frange.end() - truthHitDt));
        addTruthHit(tctr - 0.5*truthHitDt, 2.0);
        addTruthHit(tctr + 0.5*truthHitDt, 2.0);
        DOMAINCOL truthDomains; // field-off: no domains
        PTRAJPTR truthTraj = std::make_unique<PTRAJ>(truthpt);
        KKSTRAWHITCOL nostrawhits; KKSTRAWXINGCOL nostrawxings; KKCALOHITCOL nocalohits;
        // one iteration, no convergence churn: the ParameterHits already pin the params
        auto truthConfig = config_;
        truthConfig.maxniter_ = 1;
        truthConfig.divdchisq_ = std::numeric_limits<double>::max();
        truthConfig.pdchisq_ = std::numeric_limits<double>::max();
        truthConfig.divgap_ = std::numeric_limits<double>::max();
        // Extrapolation must be IDENTICAL to the reco track's. The reco fit is extended with the
        // ExtensionSettings config before extrapolating, so its extrapolation field handling is governed by
        // exconfig_.bfcorr_ (resample BField vs freeze BNom) and exconfig_.tol_ (BField-correction tolerance).
        // Carry those two fields onto the truth track's config so the extrapolation machinery is identical --
        // WITHOUT re-running the full extension schedule, which diverges a hit-less truth-pinned fit. (KKLine
        // is field-off so this is a no-op here, but keeps the two fitters consistent.)
        truthConfig.bfcorr_ = exconfig_.bfcorr_;
        truthConfig.tol_    = exconfig_.tol_;
        auto ktrk_truth = std::make_unique<KKTRK>(truthConfig,*kkbf_,kseed.particle(),truthTraj,
            nostrawhits,nostrawxings,nocalohits,truthParamHits,truthDomains);
        if(ktrk_truth->fitStatus().usable()){
          if(extrap_) extrap_->extrapolate(*ktrk_truth);
          // Sample the (extrapolated) truth track at the configured surfaces (tracker ends + CRV sectors),
          // exactly as the reco does (via RegrowKalSeed's kkfit_.sampleFit). The truth track has no straw
          // hits, so without this the CRV crossings are not recorded -- leaving trktruthsegs empty at the CRV.
          kkfit_.sampleFit(*ktrk_truth);
          TrkFitFlag tflag(kseed.status()); tflag.merge(TrkFitFlag::KKLine);
          auto tsseed = kkfit_.createSeed(*ktrk_truth,tflag,*calo_h,*nominalTracker_h);
          // carry the reco fit's hit metadata onto the truth KalSeed (the truth-constrained fit has no
          // detector hits, but EventNtuple consumers use the hit list for bookkeeping)
          tsseed._hits = kseed.hits();
          tsseed._hitcalibs = kseed.hitCalibInfos();
          tsseed._chit = kseed.caloHit();
          tsseedcol->push_back(std::move(tsseed));
          // associate the truth-seeded KalSeed to the same KalSeedMC as its reco parent
          auto tsseedp = art::Ptr<KalSeed>(KalSeedCollectionPID,tsseedcol->size()-1,KalSeedCollectionGetter);
          tsmca->addSingle(tsseedp,mcseedp);
          ktrkcol->push_back(ktrk_truth.release());
        } else if(debug_ > 0){
          std::cout << "KinematicLineTruthSeed ParameterHit fit unusable: " << ktrk_truth->fitStatus() << std::endl;
        }
      } catch (std::exception const& error) {
        if(debug_ > 0) std::cout << "KinematicLineTruthSeed error: " << error.what() << std::endl;
      }
    }
    if(debug_ > 0) std::cout << "KinematicLineTruthSeed made " << tsseedcol->size()
      << " truth-seeded tracks from " << kseedcol.size() << " reco KalSeeds" << std::endl;
    event.put(std::move(ktrkcol));
    event.put(std::move(tsseedcol));
    event.put(std::move(tsmca));
  }
}
DEFINE_ART_MODULE(mu2e::KinematicLineTruthSeed)
