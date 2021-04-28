//
// KinKal fit module using the KinematicLine parameterset
//
// Original author S. Middleton (Caltech) 2021
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
#include "RecoDataProducts/inc/CaloCluster.hh"
// data
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/KKLine.hh" 
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/Trajectory/Line.hh"
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

  using KTRAJ= KinKal::KinematicLine;
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

  class NewKinematicLineFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
     
    struct ModuleSettings {
      fhicl::Sequence<art::InputTag> cosmicTrackSeedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Cosmic seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Sequence<std::string> lineFlags { Name("LineFlags"), Comment("Flags required to be present to convert a cosmic track seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Atom<bool> saveFull { Name("SaveFullFit"), Comment("Save all track segments associated with the fit"), false};
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices"), std::vector<float>()};
      fhicl::Sequence<std::string> addHitFlags { Name("AddHitFlags"), Comment("Flags required to be present to add a hit"), std::vector<std::string>() };
      fhicl::Sequence<std::string> rejectHitFlags { Name("RejectHitFlags"), Comment("Flags required not to be present to add a hit"), std::vector<std::string>() };
      fhicl::Atom<float> maxAddDOCA { Name("MaxAddDOCA"), Comment("Max DOCA to add a hit (mm)"), 2.75 };
      fhicl::Atom<float> maxAddDt { Name("MaxAddDt"), Comment("Max Detla time to add a hit (ns)"), 3.0 };
      fhicl::Atom<float> maxAddChi { Name("MaxAddChi"), Comment("Max Chi to add a hit"), 4.0 };
      fhicl::Atom<float> maxAddDeltaU { Name("MaxAddDeltaU"), Comment("Max Delta-U to add a hit (mm)"), 10.0 };
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
    explicit NewKinematicLineFit(const ModuleParams& config);
    virtual ~NewKinematicLineFit();
    void beginRun(art::Run& run) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    KTRAJ makeSeedTraj(CosmicTrackSeed const& hseed) const;
    // data payload
    std::vector<art::ProductToken<CosmicTrackSeedCollection>> hseedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<StrawHitFlagCollection> shfcol_T_;
    TrkFitFlag goodline_;
    bool saveall_, savefull_;
    std::vector<float> zsave_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    int print_;
    float maxDoca_, maxDt_, maxChi_, maxDU_;
    KKFIT kkfit_; // fit helper	
    KKMaterial kkmat_; // material helper
    DMAT seedcov_; // seed covariance matrix
    double mass_; // particle mass
    int charge_; // particle charge
    std::unique_ptr<KKBField> kkbf_;
  };

  NewKinematicLineFit::NewKinematicLineFit(const ModuleParams& config) : art::EDProducer{config}, 
    chcol_T_(consumes<ComboHitCollection>(config().modSettings().comboHitCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(config().modSettings().strawHitFlagCollection())),
    goodline_(config().modSettings().lineFlags()),
    saveall_(config().modSettings().saveAll()),
    savefull_(config().modSettings().saveFull()),
    zsave_(config().modSettings().zsave()),
    print_(config().modSettings().printLevel()),
    maxDoca_(config().modSettings().maxAddDOCA()),
    maxDt_(config().modSettings().maxAddDt()),
    maxChi_(config().modSettings().maxAddChi()),
    maxDU_(config().modSettings().maxAddDeltaU()),
    kkfit_(config().fitSettings()),
    kkmat_(config().matSettings())
  {
    // test: only 1 of saveFull and zsave should be set
    if((savefull_ && zsave_.size() > 0) || ((!savefull_) && zsave_.size() == 0))
      throw cet::exception("RECO")<<"mu2e::NewKinematicLineFit:Segment saving configuration error"<< endl;
    // collection handling
    for(const auto& hseedtag : config().modSettings().cosmicTrackSeedCollections()) { hseedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(hseedtag)); }
    produces<KKLineCollection>();
    produces<KalSeedCollection>();
    // construct the fit configuration object.  This controls all the global and iteration-specific aspects of the fit
    // build the initial seed covariance
    auto const& seederrors = config().fitSettings().seederrors();
    if(seederrors.size() != KinKal::NParams()) 
      throw cet::exception("RECO")<<"mu2e::NewKinematicLineFit:Seed error configuration error"<< endl;
    for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
      seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];//TODO - how to turn covarience into this
    }
    // set the hit updating
    auto& schedule = kkfit_.config().schedule();
    for(auto const& shusetting : config().shuconfig() ) {
      KKStrawHitUpdater shupdater(std::get<0>(shusetting), std::get<1>(shusetting), kkfit_.nullDimension());
      unsigned minmeta = std::get<2>(shusetting);
      unsigned maxmeta = std::get<3>(shusetting);
      if(maxmeta < minmeta || schedule.size() < maxmeta)
	throw cet::exception("RECO")<<"mu2e::NewKinematicLineFit: Hit updater configuration error"<< endl;
      for(unsigned imeta=minmeta; imeta<=maxmeta; imeta++)
	schedule[imeta].updaters_.push_back(shupdater);
    }
    if(print_ > 0) std::cout << kkfit_.config();
  }

  NewKinematicLineFit::~NewKinematicLineFit(){}

  void NewKinematicLineFit::beginRun(art::Run& run) {
    // setup particle parameters
    auto const& ptable = GlobalConstantsHandle<ParticleDataTable>();
    mass_ = ptable->particle(kkfit_.fitParticle()).ref().mass().value(); 
    charge_ = static_cast<int>(ptable->particle(kkfit_.fitParticle()).ref().charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void NewKinematicLineFit::produce(art::Event& event ) {
    GeomHandle<mu2e::Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKLineCollection> kktrkcol(new KKLineCollection );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection ); //Needs to return a KalSeed
    // find the track seed collections
    for (auto const& hseedtag : hseedCols_) {
      auto const& hseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];
        
        art::Ptr<HelixSeed> hptr;

        // check helicity.  The test on the charge and helicity 
        if(hseed.status().hasAllProperties(goodline_) ){
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
          if (kkfit_.useCalo() && hptr->caloCluster())kkfit_.makeCaloHit(hptr->caloCluster(),*calo_h, pseedtraj, hits);
          if(print_ > 0){
            std::cout << hits.size() << " Hits and " << exings.size() << " Material Xings in fit" << std::endl;
          }
          if(print_ > 2){
            for(auto const& thit : hits) thit->print(std::cout,2);
            for(auto const& exing :exings) exing->print(std::cout,1);
          }
          // set the seed range given the hit TPOCA values
          seedtraj.range() = kkfit_.range(hits,exings);
          if(print_ > 0){
            //std::cout << "Seed line parameters " << hseed.track() << std::endl;
            seedtraj.print(std::cout,print_);
          }
          // create and fit the track  
          auto kktrk = make_unique<KKTRK>(kkfit_.config(),*kkbf_,seedtraj,hits,exings); //TODO - check on this
          bool save(false);
          if(save || saveall_){
            // convert fits into CosmicKalSeeds for persistence	
            kkseedcol->push_back(kkfit_.createSeed(*kktrk,hptr,zsave_,savefull_)); //TODO
            kkseedcol->back()._status.merge(TrkFitFlag::KKLine);
            kktrkcol->push_back(kktrk.release());
          }
        }
      }
    }
    // put the output products into the event
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
  }

  KTRAJ NewKinematicLineFit::makeSeedTraj(CosmicTrackSeed const& hseed) const {
    //exctract CosmicTrack (contains parameters)
    auto const& scosmic = hseed.track();
    VEC3 bnom(0.0,0.0,0.0);
    // create a PKTRAJ from the CosmicTrack fit result, to seed the KinKal fit.  First, translate the parameters
    std::tuple <double, double, double, double, double, double> info = scosmic.KinKalTrackParams();//d0,phi0,z0,cost,t0,mom
    DVEC pars;
    pars[KTRAJ::d0_] = get<0>(info);
    pars[KTRAJ::phi0_] = get<1>(info);
    pars[KTRAJ::z0_] = get<2>(info);
    pars[KTRAJ::cost_] = get<3>(info);
    pars[KTRAJ::t0_] = get<4>(info); //TODO
    pars[KTRAJ::mom_] = get<5>(info); //TODO

    // create the initial trajectory
    Parameters kkpars(pars,seedcov_); //TODO seedcov
    //  construct the seed trajectory 
    return KTRAJ(kkpars, mass_, charge_, bnom, TimeRange()); //TODO: better constructor
  }
}
DEFINE_ART_MODULE(mu2e::NewKinematicLineFit);
