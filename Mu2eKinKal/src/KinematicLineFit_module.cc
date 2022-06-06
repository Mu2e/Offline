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
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KKLine.hh"
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/Trajectory/Line.hh"
// Mu2eKinKal
#include "Offline/Mu2eKinKal/inc/KKFit.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
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
  using KTRAJ= KinKal::KinematicLine;
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
  using HPtr = art::Ptr<CosmicTrackSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using Mu2eConfig = Mu2eKinKal::Mu2eConfig;
  using KKMaterialConfig = KKMaterial::Config;

  class KinematicLineFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    struct ModuleConfig {
      fhicl::Sequence<art::InputTag> cosmicTrackSeedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Cosmic seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Sequence<std::string> lineFlags { Name("LineFlags"), Comment("Flags required to be present to convert a cosmic track seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Atom<bool> saveFull { Name("SaveFullFit"), Comment("Save all track segments associated with the fit"), false};
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices"), std::vector<float>()};
      fhicl::Sequence<std::string> addHitFlags { Name("AddHitFlags"), Comment("Flags required to be present to add a hit"), std::vector<std::string>() };
    };

    struct GlobalConfig {
      fhicl::Table<ModuleConfig> modSettings { Name("ModuleSettings") };
      fhicl::Table<Mu2eConfig> mu2eSettings { Name("Mu2eSettings") };
      fhicl::Table<KKConfig> kkFitSettings { Name("KinKalFitSettings") };
      fhicl::Table<KKConfig> kkExtSettings { Name("KinKalExtensionSettings") };
      fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
      //      StrawHitUpdateSettings shuconfig { Name("StrawHitUpdateSettings"), Comment("Setting sequence for updating StrawHits, format: \n"
      //      " 'MinDoca', 'MaxDoca', First Meta-iteration', 'Last Meta-iteration'") };
    };
    using GlobalSettings = art::EDProducer::Table<GlobalConfig>;

    public:
    explicit KinematicLineFit(const GlobalSettings& settings);
    virtual ~KinematicLineFit();
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
    Config config_; // initial fit configuration object
    Config exconfig_; // extension configuration object
  };

  KinematicLineFit::KinematicLineFit(const GlobalSettings& settings) : art::EDProducer{settings},
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(settings().modSettings().strawHitFlagCollection())),
    goodline_(settings().modSettings().lineFlags()),
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
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit:Segment saving configuration error"<< endl;
      // collection handling
      for(const auto& hseedtag : settings().modSettings().cosmicTrackSeedCollections()) { hseedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(hseedtag)); }
      produces<KKLineCollection>();
      produces<KalSeedCollection>();
      produces<KalLineAssns>();
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      if(print_ > 0) std::cout << config_;


    }

  KinematicLineFit::~KinematicLineFit(){}

  void KinematicLineFit::beginRun(art::Run& run) {
    // setup particle parameters
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(kkfit_.fitParticle()).mass();
    charge_ = static_cast<int>(ptable->particle(kkfit_.fitParticle()).charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void KinematicLineFit::produce(art::Event& event ) {
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
    unique_ptr<KalLineAssns> kkseedassns(new KalLineAssns());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // find the track seed collections
    for (auto const& hseedtag : hseedCols_) {
      auto const& hseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];

        auto hptr = HPtr(hseedcol_h,iseed);
        // check helicity.  The test on the charge and helicity
        if(hseed.status().hasAllProperties(goodline_) ){
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

          //here
          KKCALOHITCOL calohits;
          //if (kkfit_.useCalo()) kkfit_.makeCaloHit(hptr->caloCluster(),*calo_h, pseedtraj, calohits); --> CosmicTrackSeed has no CaloClusters....

          if(print_ > 2){
            for(auto const& strawhit : strawhits) strawhit->print(std::cout,2);
            for(auto const& calohit : calohits) calohit->print(std::cout,2);
            for(auto const& strawxing :strawxings) strawxing->print(std::cout,2);
          }
          // set the seed range given the hit TPOCA values
          seedtraj.range() = kkfit_.range(strawhits,calohits, strawxings);
          if(print_ > 0){
            //std::cout << "Seed line parameters " << hseed.track() << std::endl;
            seedtraj.print(std::cout,print_);
          }
          // create and fit the track
          auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,kkfit_.fitParticle(),strawhits,calohits,strawxings); //TODO - check on this
          bool save(true);//TODO - when would we like not to save?
          if(save || saveall_){
            // convert KKTrk into KalSeeds for persistence
            auto const& fittraj = kktrk->fitTraj();
            // convert fits into CosmicKalSeeds for persistence
            TrkFitFlag fitflag(hptr->status());
            fitflag.merge(TrkFitFlag::KKLine);
            // Decide which segments to save
            std::set<double> savetimes;
            if(savefull_){
              // loop over all pieces of the fit trajectory and record their times
              for (auto const& traj : fittraj.pieces() ) savetimes.insert(traj->range().mid());
            } else {
              for(auto zpos : zsave_ ) {
                // compute the time the trajectory crosses this plane
                double tz = kkfit_.zTime(fittraj,zpos);
                // find the explicit trajectory piece at this time, and store the midpoint time.  This enforces uniqueness (no duplicates)
                auto const& zpiece = fittraj.nearestPiece(tz);
                savetimes.insert(zpiece.range().mid());
              }
            }

            kkseedcol->push_back(kkfit_.createSeed(*kktrk,fitflag,*calo_h,savetimes));
            //kkseedcol->back()._status.merge(TrkFitFlag::KKLine);
            kktrkcol->push_back(kktrk.release());
            // fill assns with the cosmic seed
            auto hptr = art::Ptr<CosmicTrackSeed>(hseedcol_h,iseed);
            auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
            kkseedassns->addSingle(kseedptr,hptr);
            // save (unpersistable) KKTrk in the event
          }
        }
      }
    }
    // put the output products into the event
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
  }

  KTRAJ KinematicLineFit::makeSeedTraj(CosmicTrackSeed const& hseed) const {
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
DEFINE_ART_MODULE(mu2e::KinematicLineFit);
