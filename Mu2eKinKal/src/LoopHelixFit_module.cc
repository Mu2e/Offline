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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
// conditions
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "HepPDT/ParticleData.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
// utiliites
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeneralUtilities/inc/Angles.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
// data
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/KKLoopHelix.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/General/Parameters.hh"
// Mu2eKinKal
#include "Mu2eKinKal/inc/KKFit.hh"
#include "Mu2eKinKal/inc/KKTrack.hh"
#include "Mu2eKinKal/inc/KKMaterial.hh"
#include "Mu2eKinKal/inc/KKStrawHit.hh"
#include "Mu2eKinKal/inc/KKStrawXing.hh"
#include "Mu2eKinKal/inc/KKCaloHit.hh"
#include "Mu2eKinKal/inc/KKBField.hh"
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
  using KKFIT = KKFit<KTRAJ>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKFIT = mu2e::KKFit<KTRAJ>;
  using KinKal::DVEC;
  using KinKal::Parameters;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::DMAT;
  using KinKal::Status;
  using HPtr = art::Ptr<HelixSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;
  
  using KKConfig = Mu2eKinKal::KKConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKMaterialConfig = KKMaterial::Config;

  class LoopHelixFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

// direct configuration used explicitly in this module.  Other configurations are handled in subclasses
    struct ModuleConfig {
      fhicl::Sequence<art::InputTag> helixSeedCollections         {Name("HelixSeedCollections"),     Comment("Helix seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Sequence<std::string> helixFlags { Name("HelixFlags"), Comment("Flags required to be present to convert a helix seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };
      fhicl::Atom<bool> addHits { Name("addHits"), Comment("Add hits after initial fit"),false };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Atom<bool> saveFull { Name("SaveFullFit"), Comment("Save all helix segments associated with the fit"), false};
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices"), std::vector<float>()};
    };

    struct GlobalConfig {
      fhicl::Table<ModuleConfig> modSettings { Name("ModuleSettings") };
      fhicl::Table<KKFitConfig> kkFitSettings { Name("KKFitSettings") };
      fhicl::Table<KKConfig> kkSettings { Name("KinKalSettings") };
      fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
//      StrawHitUpdateSettings shuconfig { Name("StrawHitUpdateSettings"), Comment("Setting sequence for updating StrawHits, format: \n"
//      " 'MinDoca', 'MaxDoca', First Meta-iteration', 'Last Meta-iteration'") };
    };
    using GlobalSettings = art::EDProducer::Table<GlobalConfig>;

    public:
    explicit LoopHelixFit(const GlobalSettings& settings);
    virtual ~LoopHelixFit();
    void beginRun(art::Run& run) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    KTRAJ makeSeedTraj(HelixSeed const& hseed) const;
    // data payload
    std::vector<art::ProductToken<HelixSeedCollection>> hseedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<StrawHitFlagCollection> shfcol_T_;
    TrkFitFlag goodhelix_;
    bool addhits_, saveall_, savefull_;
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
    shfcol_T_(mayConsume<StrawHitFlagCollection>(settings().modSettings().strawHitFlagCollection())),
    goodhelix_(settings().modSettings().helixFlags()),
    addhits_(settings().modSettings().addHits()),
    saveall_(settings().modSettings().saveAll()),
    savefull_(settings().modSettings().saveFull()),
    zsave_(settings().modSettings().zsave()),
    print_(settings().modSettings().printLevel()),
    kkfit_(settings().kkFitSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().kkSettings()))
  {
    // test: only 1 of saveFull and zsave should be set
    if((savefull_ && zsave_.size() > 0) || ((!savefull_) && zsave_.size() == 0))
      throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Segment saving configuration error"<< endl;
    // collection handling
    for(const auto& hseedtag : settings().modSettings().helixSeedCollections()) { hseedCols_.emplace_back(consumes<HelixSeedCollection>(hseedtag)); }
    produces<KKLoopHelixCollection>();
    produces<KalSeedCollection>();
    // build the initial seed covariance
    auto const& seederrors = settings().modSettings().seederrors();
    if(seederrors.size() != KinKal::NParams()) 
      throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Seed error configuration error"<< endl;
    for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
      seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
    }
    if(print_ > 0) std::cout << config_;
  }

  LoopHelixFit::~LoopHelixFit(){}

  void LoopHelixFit::beginRun(art::Run& run) {
    // setup things that rely on data related to beginRun
    auto const& ptable = GlobalConstantsHandle<ParticleDataTable>();
    mass_ = ptable->particle(kkfit_.fitParticle()).ref().mass().value(); 
    charge_ = static_cast<int>(ptable->particle(kkfit_.fitParticle()).ref().charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void LoopHelixFit::produce(art::Event& event ) {
    //
    GeomHandle<mu2e::Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto const& chcol = *ch_H;
    // find calo clusters TODO
    // create output
    unique_ptr<KKLoopHelixCollection> kktrkcol(new KKLoopHelixCollection );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
    // find the helix seed collections
    for (auto const& hseedtag : hseedCols_) {
      auto const& hseedcol_h = event.getValidHandle<HelixSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
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
	  // construct hit amd material objects from Helix hits
	  KKSTRAWHITCOL strawhits; 
	  KKCALOHITCOL calohits;
	  KKSTRAWXINGCOL strawxings;
	  // first, we need to unwind the combohits.  We use this also to find the time range
	  StrawHitIndexCollection strawHitIdxs;
	  auto const& hhits = hseed.hits();
	  for(size_t ihit = 0; ihit < hhits.size(); ++ihit ){ hhits.fillStrawHitIndices(event,ihit,strawHitIdxs); }
	  // next, build straw hits from these
	  kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, strawhits, strawxings);
	  // optionally (and if present) add the CaloCluster hit
	  // verify the cluster looks physically reasonable before adding it TODO!
	  if (kkfit_.useCalo() && hseed.caloCluster())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, calohits);
	  // 
	  if(print_ > 0)
	    std::cout << strawhits.size() << " StrawHits and " << calohits.size() << " CaloHits and " << strawxings.size() << " Straw Xings in fit" << std::endl;
	  if(print_ > 1){
	    for(auto const& strawhit : strawhits) strawhit->print(std::cout,2);
	    for(auto const& calohit : calohits) calohit->print(std::cout,2);
	    for(auto const& strawxing :strawxings) strawxing->print(std::cout,2);
	  }
	  // set the seed range given the hit TPOCA values
	  seedtraj.range() = kkfit_.range(strawhits,calohits,strawxings);
	  if(print_ > 0){
	    std::cout << "Seed Helix parameters " << hseed.helix() << std::endl;
	    seedtraj.print(std::cout,print_);
	  }
	  // create and fit the track  
	  auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,strawhits,calohits,strawxings);
	  bool save(false);
	  if(kktrk->fitStatus().usable()){
	    // Check fit for physical consistency; fit can succeed but the result can have the wrong charge
	    auto const& midtraj = kktrk->fitTraj().nearestPiece(kktrk->fitTraj().range().mid());
	    save = midtraj.Q()*midtraj.rad() > 0;
	    if(save && addhits_) {
	      KKSTRAWHITCOL addhits;
	      KKSTRAWXINGCOL addexings;
	      kkfit_.addStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), *kktrk, chcol, strawHitIdxs, addhits, addexings ); 
	      if(print_ > 0)std::cout << "Found " << addhits.size() << " StrawHits and " << addexings.size() << " Straw Xings to add " << std::endl;
	      if(print_ > 1) {
		for(auto const& strawhit : strawhits) strawhit->print(std::cout,2);
		for(auto const& strawxing :strawxings) strawxing->print(std::cout,2);
	      }
	    }
	  }
	  if(save || saveall_){
	  // convert fits into KalSeeds for persistence	
	    kkseedcol->push_back(kkfit_.createSeed(*kktrk,hptr,zsave_,savefull_));
	    kkseedcol->back()._status.merge(TrkFitFlag::KKLoopHelix);
	    kktrkcol->push_back(kktrk.release());
	  }
	}
      }
    }
    // put the output products into the event
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
  }

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed) const {
    // compute the magnetic field at the helix center.  We only want the z compontent, as the helix fit assumes B points along Z
    auto const& shelix = hseed.helix();
    double zmin = std::numeric_limits<float>::max();
    double zmax = std::numeric_limits<float>::min();
    auto const& hits = hseed.hits();
    for( auto const& hit : hits) {
      double zpos = hit.pos().z();
      zmin = std::min(zmin,zpos);
      zmax = std::max(zmax,zpos);
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
    return KTRAJ(kkpars, mass_, charge_, bnom, TimeRange());
  }

}
// mu2e

DEFINE_ART_MODULE(mu2e::LoopHelixFit);
