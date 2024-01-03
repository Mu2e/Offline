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
#include "Offline/Mu2eUtilities/inc/CosmicTrackUtils.hh"
// data
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
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
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKMaterial.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
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
  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using KKFIT = KKFit<KTRAJ>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKFIT = mu2e::KKFit<KTRAJ>;
  using KinKal::DVEC;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::DMAT;
  using KinKal::Status;
  using HPtr = art::Ptr<CosmicTrackSeed>;
  using CCPtr = art::Ptr<CaloCluster>;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;

  using KKConfig = Mu2eKinKal::KinKalConfig;
  using KKFitConfig = Mu2eKinKal::KKFitConfig;
  using KKModuleConfig = Mu2eKinKal::KKModuleConfig;
  using KKMaterialConfig = KKMaterial::Config;

  class KinematicLineFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;
  // extend the generic module configuration as needed
    struct KKLineModuleConfig : public KKModuleConfig {
      fhicl::Sequence<art::InputTag> seedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Seed fit collections to be processed ") };
      fhicl::Atom<float> seedmom { Name("SeedMomentum"), Comment("Initial momentum value")};
      fhicl::Sequence<float> paramconstraints { Name("ParameterConstraints"), Comment("Sigma of direct gaussian constraints on each parameter (0=no constraint)")};
    };

    struct GlobalConfig {
      fhicl::Table<KKLineModuleConfig> modSettings { Name("ModuleSettings") };
      fhicl::Table<KKFitConfig> mu2eSettings { Name("KKFitSettings") };
      fhicl::Table<KKConfig> fitSettings { Name("FitSettings") };
      fhicl::Table<KKConfig> extSettings { Name("ExtensionSettings") };
      fhicl::Table<KKMaterialConfig> matSettings { Name("MaterialSettings") };
    };

    public:
    using Parameters = art::EDProducer::Table<GlobalConfig>;
    explicit KinematicLineFit(const Parameters& settings);
    virtual ~KinematicLineFit();
    void beginRun(art::Run& run) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    KTRAJ makeSeedTraj(CosmicTrackSeed const& hseed) const;
    bool goodFit(KKTRK const& ktrk) const;
    // data payload
    std::vector<art::ProductToken<CosmicTrackSeedCollection>> seedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<CaloClusterCollection> cccol_T_;
    TrkFitFlag goodline_;
    bool saveall_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    int print_;
    float seedmom_;
    PDGCode::type fpart_;
    KKFIT kkfit_; // fit helper
    KKMaterial kkmat_; // material helper
    DMAT seedcov_; // seed covariance matrix
    std::array<double,KinKal::NParams()> paramconstraints_;
    double mass_; // particle mass
    int charge_; // particle charge
    std::unique_ptr<KKBField> kkbf_;
    Config config_; // initial fit configuration object
    Config exconfig_; // extension configuration object
  };

  KinematicLineFit::KinematicLineFit(const Parameters& settings) : art::EDProducer{settings},
    chcol_T_(consumes<ComboHitCollection>(settings().modSettings().comboHitCollection())),
    cccol_T_(mayConsume<CaloClusterCollection>(settings().modSettings().caloClusterCollection())),
    goodline_(settings().modSettings().seedFlags()),
    saveall_(settings().modSettings().saveAll()),
    print_(settings().modSettings().printLevel()),
    seedmom_(settings().modSettings().seedmom()),
    fpart_(static_cast<PDGCode::type>(settings().modSettings().fitParticle())),
    kkfit_(settings().mu2eSettings()),
    kkmat_(settings().matSettings()),
    config_(Mu2eKinKal::makeConfig(settings().fitSettings())),
    exconfig_(Mu2eKinKal::makeConfig(settings().extSettings()))
    {
      // collection handling
      for(const auto& seedtag : settings().modSettings().seedCollections()) { seedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(seedtag)); }
      produces<KKTRKCOL>();
      produces<KalSeedCollection>();
      produces<KalLineAssns>();
      // build the initial seed covariance
      auto const& seederrors = settings().modSettings().seederrors();
      if(seederrors.size() != KinKal::NParams())
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit:Seed error configuration error"<< endl;
      for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
        seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
      }
      auto const& tempconstraints = settings().modSettings().paramconstraints();
      if (tempconstraints.size() == 0){
        for (size_t ipar=0;ipar<KinKal::NParams();ipar++)
          paramconstraints_[ipar] = 0.0;
      }else if (tempconstraints.size() == KinKal::NParams()){
        for (size_t ipar=0;ipar<KinKal::NParams();ipar++)
          paramconstraints_[ipar] = tempconstraints[ipar];
      }else{
        throw cet::exception("RECO")<<"mu2e::KinematicLineFit: Parameter constraint configuration error"<< endl;
      }
      if(print_ > 0) std::cout << config_;


    }

  KinematicLineFit::~KinematicLineFit(){}

  void KinematicLineFit::beginRun(art::Run& run) {
    // setup particle parameters
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
    mass_ = ptable->particle(fpart_).mass();
    charge_ = static_cast<int>(ptable->particle(fpart_).charge());
    // create KKBField
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::make_unique<KKBField>(*bfmgr,*det);
  }

  void KinematicLineFit::produce(art::Event& event ) {
    GeomHandle<mu2e::Calorimeter> calo_h;
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto cc_H = event.getValidHandle<CaloClusterCollection>(cccol_T_);
    auto const& chcol = *ch_H;
    // create output
    unique_ptr<KKTRKCOL> kktrkcol(new KKTRKCOL );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection ); //Needs to return a KalSeed
    unique_ptr<KalLineAssns> kkseedassns(new KalLineAssns());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // find the track seed collections
    unsigned nseed(0);
    for (auto const& hseedtag : seedCols_) {
      auto const& hseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(hseedtag);
      auto const& hseedcol = *hseedcol_h;
      nseed += hseedcol.size();
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
          auto chcolptr = hseed.hits().fillStrawHitIndices(strawHitIdxs, StrawIdMask::uniquestraw);
//          if(chcolptr != &chcol)
//            throw cet::exception("RECO")<<"mu2e::KinematicLineFit: inconsistent ComboHitCollection" << std::endl;
          // next, build straw hits and materials from these
          KKSTRAWHITCOL strawhits;
          KKSTRAWXINGCOL strawxings;
          strawhits.reserve(strawHitIdxs.size());
          strawxings.reserve(strawHitIdxs.size());
          kkfit_.makeStrawHits(*tracker, *strawresponse, *kkbf_, kkmat_.strawMaterial(), pseedtraj, *chcolptr, strawHitIdxs, strawhits, strawxings);

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
          auto kktrk = make_unique<KKTRK>(config_,*kkbf_,seedtraj,fpart_,kkfit_.strawHitClusterer(),strawhits,strawxings,calohits,paramconstraints_);
          auto goodfit = goodFit(*kktrk);
          if(goodfit && exconfig_.schedule().size() > 0){
            kkfit_.extendTrack(exconfig_,*kkbf_, *tracker,*strawresponse, kkmat_.strawMaterial(), chcol, *calo_h, cc_H, *kktrk );
          }
          bool save = goodFit(*kktrk);
          if(save || saveall_){
            TrkFitFlag fitflag(hptr->status());
            fitflag.merge(TrkFitFlag::KKLine);
            kkseedcol->push_back(kkfit_.createSeed(*kktrk,fitflag,*calo_h));
            kkseedcol->back()._status.merge(TrkFitFlag::KKLine);
            // fill assns with the cosmic seed
            auto hptr = art::Ptr<CosmicTrackSeed>(hseedcol_h,iseed);
            auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kkseedcol->size()-1,KalSeedCollectionGetter);
            kkseedassns->addSingle(kseedptr,hptr);
            // save (unpersistable) KKTrk in the event
            kktrkcol->push_back(kktrk.release());
          }
        }
      }
    }
    // put the output products into the event
    if(print_ > 0) std::cout << "Fitted " << kktrkcol->size() << " tracks from " << nseed << " seeds" << std::endl;
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
    event.put(move(kkseedassns));
  }

  KTRAJ KinematicLineFit::makeSeedTraj(CosmicTrackSeed const& hseed) const {
    //exctract CosmicTrack (contains parameters)
    VEC3 bnom(0.0,0.0,0.0);
    KinKal::VEC4 pos(hseed._track.MinuitParams.A0, 0, hseed._track.MinuitParams.B0, hseed._t0._t0);
    XYZVectorF mom3(hseed._track.MinuitParams.A1, -1, hseed._track.MinuitParams.B1);
    mom3 = mom3.Unit()*seedmom_;
    KinKal::MOM4 mom(mom3.x(),mom3.y(),mom3.z(),mass_);

    auto seedtraj = KTRAJ(pos,mom,charge_,bnom,Mu2eKinKal::timeBounds(hseed.hits()));
    seedtraj.params().covariance() = seedcov_;
    return seedtraj;
  }
  bool KinematicLineFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    return ktrk.fitStatus().usable();
  }

}
DEFINE_ART_MODULE(mu2e::KinematicLineFit)
