//
// KinKal fit module using the LoopHelix parameterset
//
// Original author D. Brown (LBNL) 11/18/20
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Tuple.h"
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
// utiliites
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/KKLoopHelix.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/General/Parameters.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
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
using namespace std;
//using namespace KinKal;
namespace mu2e {
  using KTRAJ= KinKal::LoopHelix;
  using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
  using KKTRK = KKTrack<KTRAJ>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;
  using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
  using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
  using KKSTRAWXING = KKStrawXing<KTRAJ>;
  using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
  using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
  using KKCALOHIT = KKCaloHit<KTRAJ>;
  using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
  using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using KinKal::DVEC;
  using KinKal::Parameters;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::DMAT;
  using KinKal::Status;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eConfig = Mu2eKinKal::Mu2eConfig;
  using KKMaterialConfig = KKMaterial::Config;
  using Name    = fhicl::Name;
  using Comment = fhicl::Comment;
  // direct configuration used explicitly in this module.  Other configurations are handled in subclasses
  struct ModuleConfig {
    fhicl::Sequence<art::InputTag> helixSeedCollections         {Name("HelixSeedCollections"),     Comment("Helix seed fit collections to be processed ") };
    fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
    fhicl::Atom<art::InputTag>     caloClusterCollection     {Name("CaloClusterCollection"),     Comment("CaloCluster collection ") };
    fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
    fhicl::Sequence<std::string> helixFlags { Name("HelixFlags"), Comment("Flags required to be present to convert a helix seed to a KinKal track") };
    fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
    fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };
    fhicl::Atom<bool> extend { Name("Extend"), Comment("Extend fit with hits and materials found after initial fit"),false };
    fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
    fhicl::Atom<bool> saveFull { Name("SaveFullFit"), Comment("Save all helix segments associated with the fit"), false};
    fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices"), std::vector<float>()};
  };

  struct GlobalConfig {
    fhicl::Table<ModuleConfig> modSettings { Name("ModuleSettings") };
    fhicl::Table<Mu2eConfig> mu2eSettings { Name("Mu2eSettings") };
    fhicl::Table<KKConfig> kkFitSettings { Name("KinKalFitSettings") };
    fhicl::Table<KKConfig> kkExtSettings { Name("KinKalExtensionSettings") };
    fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
  };

  class LoopHelixFit : public art::EDProducer {
    using GlobalSettings = art::EDProducer::Table<GlobalConfig>;

    public:
    explicit LoopHelixFit(const GlobalSettings& settings);
    virtual ~LoopHelixFit();
    void beginRun(art::Run& run) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    KTRAJ makeSeedTraj(HelixSeed const& hseed) const;
    bool goodFit(KKTRK const& ktrk) const;
    void fillSaveTimes(KKTRK const& ktrk,std::set<double>& savetimes) const;
    // data payload
    std::vector<art::ProductToken<HelixSeedCollection>> hseedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<CaloClusterCollection> cccol_T_;
    art::ProductToken<StrawHitFlagCollection> shfcol_T_;
    TrkFitFlag goodhelix_;
    bool extend_, saveall_, savefull_;
    std::vector<float> zsave_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    int print_;
    KKFIT kkfit_; // fit helper
    KKMaterial kkmat_; // material helper
    DMAT seedcov_; // seed covariance matrix
    double mass_; // particle mass
    int charge_; // particle charge
    std::unique_ptr<KKBField> kkbf_;
    Config config_; // initial fit configuration object
    Config exconfig_; // extension configuration object
  };

  LoopHelixFit::LoopHelixFit(const GlobalSettings& settings) : art::EDProducer{settings},
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(settings().modSettings().strawHitFlagCollection())),
    goodhelix_(settings().modSettings().helixFlags()),
    extend_(settings().modSettings().extend()),
    saveall_(settings().modSettings().saveAll()),
    savefull_(settings().modSettings().saveFull()),
    zsave_(settings().modSettings().zsave()),
    print_(settings().modSettings().printLevel()),
    kkfit_(settings().mu2eSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().kkFitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().kkExtSettings()))
    {
      // test: only 1 of saveFull and zsave should be set
      if((savefull_ && zsave_.size() > 0) || ((!savefull_) && zsave_.size() == 0))
        throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Segment saving configuration error"<< endl;
      // collection handling
      for(const auto& hseedtag : settings().modSettings().helixSeedCollections()) { hseedCols_.emplace_back(consumes<HelixSeedCollection>(hseedtag)); }
      produces<KKLoopHelixCollection>();
      produces<KalSeedCollection>();
      produces<KalHelixAssns>();
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      if(print_ > 0) std::cout << "Fit " << config_ << "Extension " << exconfig_;
    }

  LoopHelixFit::~LoopHelixFit(){}

  void LoopHelixFit::beginRun(art::Run& run) {
    // setup things that rely on data related to beginRun
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(kkfit_.fitParticle()).mass();
    charge_ = static_cast<int>(ptable->particle(kkfit_.fitParticle()).charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void LoopHelixFit::produce(art::Event& event ) {
    GeomHandle<mu2e::Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getValidHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKLoopHelixCollection> kktrkcol(new KKLoopHelixCollection );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
    unique_ptr<KalHelixAssns> kkseedassns(new KalHelixAssns());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // find the helix seed collections
    unsigned nhelix(0);
    for (auto const& hseedtag : hseedCols_) {
      auto const& hseedcol_h = event.getValidHandle<HelixSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      nhelix += hseedcol.size();
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];
        auto hptr = HPtr(hseedcol_h,iseed);
        // check helicity.  The test on the charge and helicity
        if(hseed.status().hasAllProperties(goodhelix_) ){
          // construt the seed trajectory
          KTRAJ seedtraj = makeSeedTraj(hseed);
          // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
          PKTRAJ pseedtraj(seedtraj);
          // first, we need to unwind the combohits.  We use this also to find the time range
          StrawHitIndexCollection strawHitIdxs;
          auto const& hhits = hseed.hits();
          for(size_t ihit = 0; ihit < hhits.size(); ++ihit ){ hhits.fillStrawHitIndices(event,ihit,strawHitIdxs); }
          // next, build straw hits and materials from these
          KKSTRAWHITCOL strawhits;
          KKSTRAWXINGCOL strawxings;
          strawhits.reserve(hhits.size());
          strawxings.reserve(hhits.size());
          kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings);
          // optionally (and if present) add the CaloCluster hit
          // verify the cluster looks physically reasonable before adding it TODO!  Or, let the KKCaloHit updater do it
          KKCALOHITCOL calohits;
          if (kkfit_.useCalo() && hseed.caloCluster().isNonnull())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
          // extend the seed range given the hits and xings
          seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
          // create and fit the track
          auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,kkfit_.fitParticle(),strawhits,calohits,strawxings);
          // Check the fit
          auto goodfit = goodFit(*kktrk);
          if(goodfit && extend_) {
            kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *kktrk );
            goodfit = goodFit(*kktrk);
          }
          if(print_>0)kktrk->printFit(std::cout,print_);
          if(goodfit || saveall_){
            TrkFitFlag fitflag(hptr->status());
            fitflag.merge(TrkFitFlag::KKLoopHelix);
            // Decide which segments to save
            std::set<double> savetimes;
            fillSaveTimes(*kktrk,savetimes);
            kkseedcol->push_back(kkfit_.createSeed(*kktrk,fitflag,savetimes));
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
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << nhelix << " Helices" << std::endl;
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
  }

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed) const {
    // compute the magnetic field at the helix center.  We only want the z compontent, as the helix fit assumes B points along Z
    auto const& shelix = hseed.helix();
    double zmin = std::numeric_limits<float>::max();
    double zmax = std::numeric_limits<float>::min();
    double tmin = std::numeric_limits<float>::max();
    double tmax = std::numeric_limits<float>::min();
    auto const& hits = hseed.hits();
    for( auto const& hit : hits) {
      double zpos = hit.pos().z();
      zmin = std::min(zmin,zpos);
      zmax = std::max(zmax,zpos);
      tmin = std::min(tmin,(double)hit.correctedTime());
      tmax = std::max(tmax,(double)hit.correctedTime());
    }
    float zcent = 0.5*(zmin+zmax);
    VEC3 center(shelix.centerx(), shelix.centery(),zcent);
    auto bcent = kkbf_->fieldVect(center);
    VEC3 bnom(0.0,0.0,bcent.Z());
    // create a PKTRAJ from the helix fit result, to seed the KinKal fit.  First, translate the parameters
    // Note the sign adjustments; RobustHelix is a purely geometric helix, with slightly different conventions
    DVEC pars;
    double psign = copysign(1.0,-charge_*bcent.Z());
    pars[KTRAJ::rad_] = shelix.radius()*psign;
    pars[KTRAJ::lam_] = shelix.lambda()*kkfit_.fitDirection().dzdt();
    pars[KTRAJ::cx_] = shelix.centerx();
    pars[KTRAJ::cy_] = shelix.centery();
    pars[KTRAJ::phi0_] = shelix.fz0()+psign*M_PI_2;
    pars[KTRAJ::t0_] = hseed.t0().t0();
    // create the initial trajectory
    Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    return KTRAJ(kkpars, mass_, charge_, bnom, TimeRange(tmin,tmax));
  }

  bool LoopHelixFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    return ktrk.fitStatus().usable() &&
      ktrk.fitTraj().front().Q()*ktrk.seedTraj().Q() > 0 &&
      ktrk.fitTraj().front().helicity()*ktrk.seedTraj().helicity() > 0;
  }

  void LoopHelixFit::fillSaveTimes(KKTRK const& ktrk,std::set<double>& savetimes) const {
    auto const& fittraj = ktrk.fitTraj();
    if(savefull_){ // loop over all pieces of the fit trajectory and record their times
      for (auto const& traj : fittraj.pieces() ) savetimes.insert(traj.range().mid());
    } else {
      for(auto zpos : zsave_ ) {
        // compute the time the trajectory crosses this plane
        double tz = kkfit_.zTime(fittraj,zpos);
        // find the explicit trajectory piece at this time, and store the midpoint time.  This enforces uniqueness (no duplicates)
        auto const& zpiece = fittraj.nearestPiece(tz);
        savetimes.insert(zpiece.range().mid());
      }
    }
  }

} // mu2e

DEFINE_ART_MODULE(mu2e::LoopHelixFit);
