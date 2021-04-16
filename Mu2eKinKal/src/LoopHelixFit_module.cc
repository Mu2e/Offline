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
#include "Mu2eKinKal/inc/KKMaterial.hh"
#include "Mu2eKinKal/inc/KKStrawHit.hh"
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
  using KKTRK = KinKal::Track<KTRAJ>;
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
  
  using StrawHitUpdateSettings = fhicl::Sequence<fhicl::Tuple<float,float,unsigned,unsigned>>;
  using KKFitSettings = Mu2eKinKal::FitSettings;
  using KKMaterialSettings = KKMaterial::Settings;

  class LoopHelixFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    struct ModuleSettings {
      fhicl::Sequence<art::InputTag> helixSeedCollections         {Name("HelixSeedCollections"),     Comment("Helix seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Sequence<std::string> helixFlags { Name("HelixFlags"), Comment("Flags required to be present to convert a helix seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Atom<bool> addHits { Name("addHits"), Comment("Add hits after initial fit"),false };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Atom<bool> saveFull { Name("SaveFullFit"), Comment("Save all helix segments associated with the fit"), false};
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices"), std::vector<float>()};
    };

    struct ModuleConfig {
      fhicl::Table<ModuleSettings> modSettings { Name("ModuleSettings") };
      fhicl::Table<KKFitSettings> fitSettings { Name("FitSettings") };
      fhicl::Table<KKMaterialSettings> matSettings { Name("MaterialSettings") };
      StrawHitUpdateSettings shuconfig { Name("StrawHitUpdateSettings"), Comment("Setting sequence for updating StrawHits, format: \n"
      " 'MinDoca', 'MaxDoca', First Meta-iteration', 'Last Meta-iteration'") };
    };
    using ModuleParams = art::EDProducer::Table<ModuleConfig>;

    public:
    explicit LoopHelixFit(const ModuleParams& config);
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
  };

  LoopHelixFit::LoopHelixFit(const ModuleParams& config) : art::EDProducer{config}, 
    chcol_T_(consumes<ComboHitCollection>(config().modSettings().comboHitCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(config().modSettings().strawHitFlagCollection())),
    goodhelix_(config().modSettings().helixFlags()),
    addhits_(config().modSettings().addHits()),
    saveall_(config().modSettings().saveAll()),
    savefull_(config().modSettings().saveFull()),
    zsave_(config().modSettings().zsave()),
    print_(config().modSettings().printLevel()),
    kkfit_(config().fitSettings()),
    kkmat_(config().matSettings())
  {
    // test: only 1 of saveFull and zsave should be set
    if((savefull_ && zsave_.size() > 0) || ((!savefull_) && zsave_.size() == 0))
      throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Segment saving configuration error"<< endl;
    // collection handling
    for(const auto& hseedtag : config().modSettings().helixSeedCollections()) { hseedCols_.emplace_back(consumes<HelixSeedCollection>(hseedtag)); }
    produces<KKLoopHelixCollection>();
    produces<KalSeedCollection>();
    // construct the fit configuration object.  This controls all the global and iteration-specific aspects of the fit
    // build the initial seed covariance
    auto const& seederrors = config().fitSettings().seederrors();
    if(seederrors.size() != KinKal::NParams()) 
      throw cet::exception("RECO")<<"mu2e::LoopHelixFit:Seed error configuration error"<< endl;
    for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
      seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
    }
    // set the hit updating
    auto& schedule = kkfit_.config().schedule();
    for(auto const& shusetting : config().shuconfig() ) {
      KKStrawHitUpdater shupdater(std::get<0>(shusetting), std::get<1>(shusetting), kkfit_.nullDimension());
      unsigned minmeta = std::get<2>(shusetting);
      unsigned maxmeta = std::get<3>(shusetting);
      if(maxmeta < minmeta || schedule.size() < maxmeta)
	throw cet::exception("RECO")<<"mu2e::LoopHelixFit: Hit updater configuration error"<< endl;
      for(unsigned imeta=minmeta; imeta<=maxmeta; imeta++)
	schedule[imeta].updaters_.push_back(shupdater);
    }
    if(print_ > 0) std::cout << kkfit_.config();
  }

  LoopHelixFit::~LoopHelixFit(){}

  void LoopHelixFit::beginRun(art::Run& run) {
    // setup particle parameters
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
	  MEASCOL hits; // polymorphic container of hits
	  EXINGCOL exings; // polymorphic container of detector element crossings
	  // first, we need to unwind the combohits.  We use this also to find the time range
	  StrawHitIndexCollection strawHitIdxs;
	  auto const& hhits = hseed.hits();
	  for(size_t ihit = 0; ihit < hhits.size(); ++ihit ){ hhits.fillStrawHitIndices(event,ihit,strawHitIdxs); }
	  // next, build straw hits from these
	  kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, chcol, strawHitIdxs, hits, exings);
	  // optionally (and if present) add the CaloCluster hit
	  // verify the cluster looks physically reasonable before adding it TODO!
	  if (kkfit_.useCalo() && hseed.caloCluster())kkfit_.makeCaloHit(hseed.caloCluster(),*calo_h, pseedtraj, hits);
	  // 
	  if(print_ > 0)
	    std::cout << hits.size() << " Hits and " << exings.size() << " Material Xings in fit" << std::endl;
	  if(print_ > 2){
	    for(auto const& thit : hits) thit->print(std::cout,2);
	    for(auto const& exing :exings) exing->print(std::cout,1);
	  }
	  // set the seed range given the hit TPOCA values
	  seedtraj.range() = kkfit_.range(hits,exings);
	  if(print_ > 0){
	    std::cout << "Seed Helix parameters " << hseed.helix() << std::endl;
	    seedtraj.print(std::cout,print_);
	  }
	  // create and fit the track  
	  auto kktrk = make_unique<KKTRK>(kkfit_.config(),*kkbf_,seedtraj,hits,exings);
	  bool save(false);
	  if(kktrk->fitStatus().usable()){
	    // Check fit for physical consistency; fit can succeed but the result can have the wrong charge
	    auto const& midtraj = kktrk->fitTraj().nearestPiece(kktrk->fitTraj().range().mid());
	    save = midtraj.Q()*midtraj.rad() > 0;
	    if(save && addhits_) {
	      kkfit_.addStrawHits(*kktrk, chcol, strawHitIdxs, *tracker, *strawresponse); 
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
