//
// KinKal fit module using the KinematicLine parameterset
//
// Original author S. Middleton (Caltech) 2021
//

// framework
#include "fhiclcpp/ParameterSet.h"
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
// data
#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "RecoDataProducts/inc/KKLine.hh" 
// KinKal
#include "KinKal/Fit/Track.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/Trajectory/Line.hh"
// Mu2eKinKal
#include "Mu2eKinKal/inc/KKFileFinder.hh"
#include "Mu2eKinKal/inc/KKStrawHit.hh"
//#include "Mu2eKinKal/inc/KKPanelHit.hh"
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
  using KTRAJ= KinKal::KinematicLine;
  using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
//  using TSH = KKStrawHit<KTRAJ>;
//  using TPH = KKPanelHit<KTRAJ>;
  using KKTRK = KinKal::Track<KTRAJ>;
  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;
  using STRAWXING = KinKal::StrawXing<KTRAJ>;
  using EXING = KinKal::ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using KKEFF = KinKal::Effect<KTRAJ>;
  using KKHIT = KinKal::HitConstraint<KTRAJ>;
  using KKMAT = KinKal::Material<KTRAJ>;
  using KKBF = KinKal::BFieldEffect<KTRAJ>;
  using KinKal::Line;
  using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
  using TCA = KinKal::ClosestApproach<KTRAJ,Line>;
  using MatEnv::MatDBInfo;
  using KKConfig = KinKal::Config;
  using KinKal::DVEC;
  using KinKal::Parameters;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::MetaIterConfig;
  using KinKal::CAHint;
  using KinKal::StrawMaterial;
  using KinKal::WireHitState;
  using KinKal::DMAT;
  using KinKal::Status;
  using HPtr = art::Ptr<CosmicTrackSeed>;

  using StrawHitIndexCollection = std::vector<StrawHitIndex>;


  class KinematicLineFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    struct ModuleSettings {
      fhicl::Sequence<art::InputTag> cosmicSeedCollections         {Name("CosmicTrackSeedCollections"),     Comment("Helix seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Sequence<std::string> lineFlags { Name("LineFlags"), Comment("Flags required to be present to convert a helix seed to a KinKal track") };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Diagnostic printout Level"), 0 };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit result helices")};
    };

    struct PatRecSettings {
      fhicl::Sequence<std::string> addHitFlags { Name("AddHitFlags"), Comment("Flags required to be present to add a hit"), std::vector<std::string>() };
      fhicl::Sequence<std::string> rejectHitFlags { Name("RejectHitFlags"), Comment("Flags required not to be present to add a hit"), std::vector<std::string>() };
      fhicl::Atom<float> maxAddDOCA { Name("MaxAddDOCA"), Comment("Max DOCA to add a hit (mm)"), 2.75 };
      fhicl::Atom<float> maxAddDt { Name("MaxAddDt"), Comment("Max Detla time to add a hit (ns)"), 3.0 };
      fhicl::Atom<float> maxAddChi { Name("MaxAddChi"), Comment("Max Chi to add a hit"), 4.0 };
      fhicl::Atom<float> maxAddDeltaU { Name("MaxAddDeltaU"), Comment("Max Delta-U to add a hit (mm)"), 10.0 };
    };

    struct MaterialSettings {
      fhicl::Atom<std::string> isotopes { Name("isotopes"), Comment("Filename for istotopes information")};
      fhicl::Atom<std::string> elements { Name("elements"), Comment("Filename for elements information") };
      fhicl::Atom<std::string> materials { Name("materials"), Comment("Filename for materials information") };
      fhicl::Atom<std::string> strawGasMaterialName{ Name("strawGasMaterialName"), Comment("strawGasMaterialName") };
      fhicl::Atom<std::string> strawWallMaterialName{ Name("strawWallMaterialName"), Comment("strawWallMaterialName") };
      fhicl::Atom<std::string> strawWireMaterialName{ Name("strawWireMaterialName"), Comment("strawWireMaterialName") };
      fhicl::Atom<double> dahlLynchScatteringFraction{ Name("dahlLynchScatteringFraction"), Comment("dahlLynchScatteringFraction") };
    };

    struct FitSettings {
      fhicl::Atom<int> fitParticle {  Name("FitParticle"), Comment("Particle type to fit: e-, e+, mu-, ..."), PDGCode::e_minus};
      fhicl::Atom<int> fitDirection { Name("FitDirection"), Comment("Particle direction to fit, either upstream or downstream"), TrkFitDirection::downstream };
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration"), 10 };
      fhicl::Atom<float> dwt { Name("Deweight"), Comment("Deweighting factor when initializing the track end parameters"), 1.0e6 };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for final fit (ns)"), 0.2 };
      fhicl::Atom<float> btol { Name("BCorrTolerance"), Comment("Tolerance on BField correction accuracy (mm)"), 0.01 };
      fhicl::Sequence<float> seederrors { Name("SeedErrors"), Comment("Initial value of seed parameter errors (rms, various units)") };//TODO
      fhicl::Atom<int> bfieldCorr { Name("BFieldCorrection"), Comment("BField correction algorithm") };
      fhicl::Atom<int> minndof { Name("MinNDOF"), Comment("Minimum number of Degrees of Freedom to conitnue fitting"), 5  };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Internal fit print level"),0};
      fhicl::Atom<int> nullHitDimension { Name("NullHitDimension"), Comment("Null hit constrain dimension"), 2 }; 
      fhicl::Atom<float> nullHitVarianceScale { Name("NullHitVarianceScale"), Comment("Scale factor on geometric null hit variance"), 1.0 }; 
      fhicl::Atom<float> tPOCAPrec { Name("TPOCAPrecision"), Comment("Precision for TPOCA calculation (ns)"), 1e-5 };
      fhicl::Atom<bool> addMaterial { Name("AddMaterial"), Comment("Add material effects to the fit"), true }; 
  };

    using MetaIterationSettings = fhicl::Sequence<fhicl::Tuple<float,float,float>>;
    using StrawHitUpdateSettings = fhicl::Sequence<fhicl::Tuple<float,float,unsigned,unsigned>>;
    struct ModuleConfig {
      fhicl::Table<ModuleSettings> modsettings { Name("ModuleSettings") };
      fhicl::Table<PatRecSettings> patrecsettings { Name("PatRecSettings") };
      fhicl::Table<FitSettings> fitsettings { Name("FitSettings") };
      fhicl::Table<MaterialSettings> matsettings { Name("MaterialSettings") };
      MetaIterationSettings mconfig { Name("MetaIterationSettings"), Comment("MetaIteration sequence configuration parameters, format: \n"
      " 'Temperature (dimensionless)', Delta chisquared/DOF for convergence', 'Delta chisquared/DOF for divergence'") };
      StrawHitUpdateSettings shuconfig { Name("StrawHitUpdateSettings"), Comment("Setting sequence for updating StrawHits, format: \n"
      " 'MinDoca', 'MaxDoca', First Meta-iteration', 'Last Meta-iteration'") };
    };
    using ModuleParams = art::EDProducer::Table<ModuleConfig>;

    public:
    explicit KinematicLineFit(const ModuleParams& config);
    virtual ~KinematicLineFit();
    void beginRun(art::Run& run) override;
    void beginSubRun(art::SubRun& subrun) override;
    void produce(art::Event& event) override;
    private:
    // utility functions
    CosmicKalSeed createSeed(KKTRK const& ktrk,HPtr const& hptr) const; //TODO - should this be KalCosmicSeed?
    KTRAJ makeSeedTraj(CosmicTrackSeed const& hseed) const;
    void makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,KTRAJ const& ktraj,ComboHitCollection const& hhits,StrawHitIndexCollection const& shidxs,
      MEASCOL& thits,EXINGCOL& exings) const;

    double zTime(PKTRAJ const& trak, double zpos) const; // find the time the trajectory crosses the plane perp to z at the given z position
    // data payload
    std::vector<art::ProductToken<CosmicTrackSeedCollection>> hseedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<StrawHitFlagCollection> shfcol_T_;
    TrkFitFlag goodline_;
    bool saveall_;
    std::vector<float> zsave_;
    TrkFitDirection tdir_;
    PDGCode::type tpart_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    int print_;
    float maxDoca_, maxDt_, maxChi_, maxDU_, tbuff_, tpocaprec_;
    WireHitState::Dimension nulldim_;
    float nullvscale_;
    bool addmat_;
    WireHitState whstate_; // initial state for new straw hits
    KKFileFinder filefinder_;
    std::string wallmatname_, gasmatname_, wirematname_;
    std::unique_ptr<StrawMaterial> smat_; // straw material
    KKConfig kkconfig_; // KinKal fit configuration
    DMAT seedcov_; // seed covariance matrix
    MatDBInfo* matdbinfo_; // material database
    double mass_; // caches
    int charge_;
    std::unique_ptr<KKBField> kkbf_;
  };

  KinematicLineFit::KinematicLineFit(const ModuleParams& config) : art::EDProducer{config}, 
    chcol_T_(consumes<ComboHitCollection>(config().modsettings().comboHitCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(config().modsettings().strawHitFlagCollection())),
    goodline_(config().modsettings().lineFlags()),
    saveall_(config().modsettings().saveAll()),
    zsave_(config().modsettings().zsave()),
    tdir_(static_cast<TrkFitDirection::FitDirection>(config().fitsettings().fitDirection())),       tpart_(static_cast<PDGCode::type>(config().fitsettings().fitParticle())),
    print_(config().modsettings().printLevel()),
    maxDoca_(config().patrecsettings().maxAddDOCA()),
    maxDt_(config().patrecsettings().maxAddDt()),
    maxChi_(config().patrecsettings().maxAddChi()),
    maxDU_(config().patrecsettings().maxAddDeltaU()),
    tbuff_(config().fitsettings().tBuffer()),
    tpocaprec_(config().fitsettings().tPOCAPrec()),
    nulldim_(static_cast<WireHitState::Dimension>(config().fitsettings().nullHitDimension())),
    nullvscale_(config().fitsettings().nullHitVarianceScale()),
    addmat_(config().fitsettings().addMaterial()),
    filefinder_(config().matsettings().elements(),config().matsettings().isotopes(),config().matsettings().materials()),
    wallmatname_(config().matsettings().strawWallMaterialName()),
    gasmatname_(config().matsettings().strawGasMaterialName()),
    wirematname_(config().matsettings().strawWireMaterialName()),
    matdbinfo_(0)
  {
    // collection handling
    for(const auto& hseedtag : config().modsettings().cosmicSeedCollections()) { hseedCols_.emplace_back(consumes<CosmicTrackSeedCollection>(hseedtag)); }
    produces<KKLineCollection>();
    produces<CosmicKalSeedCollection>();//TODO
    // construct the fit configuration object.  This controls all the global and iteration-specific aspects of the fit
    kkconfig_.maxniter_ = config().fitsettings().maxniter();
    kkconfig_.dwt_ = config().fitsettings().dwt();
    kkconfig_.tbuff_ = config().fitsettings().tBuffer();
    kkconfig_.tol_ = config().fitsettings().btol();
    kkconfig_.minndof_ = config().fitsettings().minndof();
    kkconfig_.bfcorr_ = static_cast<KKConfig::BFCorr>(config().fitsettings().bfieldCorr());
    kkconfig_.plevel_ = static_cast<KKConfig::printLevel>(config().fitsettings().printLevel());
    // build the initial seed covariance
    auto const& seederrors = config().fitsettings().seederrors(); //TODO - seed errors?
    if(seederrors.size() != KinKal::NParams()) 
      throw cet::exception("RECO")<<"mu2e::KinematicLineFit:Seed error configuration error"<< endl;
    for(size_t ipar=0;ipar < seederrors.size(); ++ipar){
      seedcov_[ipar][ipar] = seederrors[ipar]*seederrors[ipar];
    }
    // Now set the schedule for the meta-iterations
    unsigned nmiter(0);
    for(auto const& misetting : config().mconfig()) {
      MetaIterConfig mconfig;
      mconfig.temp_ = std::get<0>(misetting);
      mconfig.tprec_ = tpocaprec_;
      mconfig.convdchisq_ = std::get<1>(misetting);
      mconfig.divdchisq_ = std::get<2>(misetting);
      mconfig.miter_ = nmiter++;
      kkconfig_.schedule_.push_back(mconfig);
    }
    auto& schedule = kkconfig_.schedule();
    // simple hit updating
    for(auto const& shusetting : config().shuconfig() ) {
      KKStrawHitUpdater shupdater(std::get<0>(shusetting), std::get<1>(shusetting), static_cast<WireHitState::Dimension>(config().fitsettings().nullHitDimension()));
      unsigned minmeta = std::get<2>(shusetting);
      unsigned maxmeta = std::get<3>(shusetting);
      if(maxmeta < minmeta || kkconfig_.schedule_.size() < maxmeta)
      throw cet::exception("RECO")<<"mu2e::KinematicLineFit: Hit updater configuration error"<< endl;
      for(unsigned imeta=minmeta; imeta<=maxmeta; imeta++){
        schedule[imeta].updaters_.push_back(shupdater);
      }
    }
  }

  KinematicLineFit::~KinematicLineFit(){}

  void KinematicLineFit::beginRun(art::Run& run) {
    // initialize material access
    GeomHandle<Tracker> tracker;

    auto const& ptable = GlobalConstantsHandle<ParticleDataTable>();
    auto const& sprop = tracker->strawProperties();
    matdbinfo_ = new MatDBInfo(filefinder_);
    smat_ = std::make_unique<StrawMaterial>(
    sprop._strawOuterRadius, sprop._strawWallThickness, sprop._wireRadius,
    matdbinfo_->findDetMaterial(wallmatname_),
    matdbinfo_->findDetMaterial(gasmatname_),
    matdbinfo_->findDetMaterial(wirematname_));
  
    // kinematic properties
    mass_ = ptable->particle(tpart_).ref().mass().value(); 
    charge_ = static_cast<int>(ptable->particle(tpart_).ref().charge()); //TODO - TPART for cosmics?
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    kkbf_ = std::move(std::make_unique<KKBField>(*bfmgr,*det));
  }

  void KinematicLineFit::beginSubRun(art::SubRun& subrun) {
  }

  void KinematicLineFit::produce(art::Event& event ) {
    // find current proditions
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    // initialize hits as null (no drift)
    auto const& sprop = tracker->strawProperties();
    double rstraw = sprop.strawInnerRadius();
    double nulldt = 0.5*rstraw/strawresponse->driftConstantSpeed(); // approximate shift in time due to ignoring drift
    double nullvar = nullvscale_*rstraw*rstraw/3.0; // scaled square RMS (distance is between 0 and r)
    whstate_ = WireHitState(WireHitState::null, nulldim_, nullvar ,nulldt);
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto const& chcol = *ch_H;
    // find calo clusters TODO
    // create output
    unique_ptr<KKLineCollection> kktrkcol(new KKLineCollection );
    unique_ptr<CosmicKalSeedCollection> kkseedcol(new CosmicKalSeedCollection );
    // find the  seed collections
    for (auto const& hseedtag : hseedCols_) { //loop over seed collections 
      auto const& hseedcol_h = event.getValidHandle<CosmicTrackSeedCollection>(hseedtag); //find cosmic seeds
      auto const& hseedcol = *hseedcol_h;
      // loop over the seeds
      for(size_t iseed=0; iseed < hseedcol.size(); ++iseed) {
        auto const& hseed = hseedcol[iseed];
        auto hptr = HPtr(hseedcol_h,iseed);
        if(hseed.status().hasAllProperties(goodline_)){
          // construt the seed trajectory
          KTRAJ seedtraj = makeSeedTraj(hseed); 
          // construct hit amd material objects from Helix hits
          MEASCOL thits; // polymorphic container of hits
          EXINGCOL exings; // polymorphic container of detector element crossings
          // first, we need to unwind the combohits.  We use this also to find the time range
          StrawHitIndexCollection strawHitIdxs;
          auto const& hhits = hseed.hits();
          //loop over hits - find indicies:
          for(size_t ihit = 0; ihit < hhits.size(); ++ihit ){ hhits.fillStrawHitIndices(event,ihit,strawHitIdxs); }
          makeStrawHits(*tracker, *strawresponse, seedtraj,chcol,strawHitIdxs,thits,exings);
          // add nearby unused straws as passive material TODO
          // create and fit the track  
          auto kktrk = make_unique<KKTRK>(kkconfig_,*kkbf_,seedtraj,thits,exings);
          if(kktrk->fitStatus().usable() || saveall_){
            // convert fits into KalSeeds for persistence	
            kkseedcol->push_back(createSeed(*kktrk,hptr));
            kktrkcol->push_back(kktrk.release());
          }
        }
      }
    }
    // put the output products into the event
    event.put(move(kktrkcol));
    event.put(move(kkseedcol));
  }

  CosmicKalSeed KinematicLineFit::createSeed(KKTRK const& kktrk,HPtr const& hptr) const { //HPTR?
    //check status:
    TrkFitFlag fflag(hptr->status());
    fflag.merge(TrkFitFlag::KKLine);
    if(kktrk.fitStatus().usable()) fflag.merge(TrkFitFlag::kalmanOK);
    if(kktrk.fitStatus().status_ == Status::converged) fflag.merge(TrkFitFlag::kalmanConverged);
    // explicit T0 is needed for backwards-compatibility; sample from the appropriate trajectory piece
    auto const& fittraj = kktrk.fitTraj();
    double tz0 = zTime(fittraj,0.0);
    auto const& t0piece = fittraj.nearestPiece(tz0);//TODO
    HitT0 t0(t0piece.paramVal(KTRAJ::t0_), sqrt(t0piece.paramVal(KTRAJ::t0_))); //TODO - ParamVar
    // create the shell for the output.  Note the (obsolete) flight length is given as t0
    CosmicKalSeed fseed(tpart_,tdir_,fflag,t0.t0()); 
    fseed._cosmicseed = hptr; //TODO shouldnt be helix
    auto const& fstatus = kktrk.fitStatus();
    fseed._chisq = fstatus.chisq_.chisq();
    fseed._fitcon = fstatus.chisq_.probability();
    // loop over individual effects
    for(auto const& eff: kktrk.effects()) {
      const KKHIT* kkhit = dynamic_cast<const KKHIT*>(eff.get());
      if(kkhit != 0){
        const KKSTRAWHIT* strawhit = dynamic_cast<const KKSTRAWHIT*>(kkhit->hit().get());

        if(strawhit != 0) {
          auto const& chit = strawhit->hit();
          StrawHitFlag hflag = chit.flag();
          if(strawhit->active())hflag.merge(StrawHitFlag::active);
          auto const& ca = strawhit->closestApproach();
          TrkStrawHitSeed seedhit(strawhit->strawHitIndex(), // drift radius and drift radius error are unfilled TODO
          HitT0(ca.particleToca(),sqrt(ca.tocaVar())),
          static_cast<float>(ca.particleToca()), static_cast<float>(ca.sensorToca()),
          static_cast<float>(-1.0), static_cast<float>(strawhit->time()),
          static_cast<float>(ca.doca()), strawhit->hitState().lrambig_,static_cast<float>(-1.0), hflag, chit);
          fseed._hits.push_back(seedhit);

          // save the segment at the CH
          auto const& zpiece = fittraj.nearestPiece(ca.particleToca());
          fseed._segments.emplace_back(zpiece,ca.particleToca()); //segements
        }
	    }
      }
      //      const KKMAT* kkmat = dynamic_cast<const KKMAT*>(eff.get());
      //TrkUtilities::fillStraws(ktrk,fseed._straws); TODO!!
    // sample the fit at the requested z positions.
    for(auto zpos : zsave_ ) { //TODO - what are our save Z ?
      // compute the time the trajectory crosses this plane
      double tz = zTime(fittraj,zpos);
      // find the explicit trajectory piece at this time
      auto const& zpiece = fittraj.nearestPiece(tz);
      // construct and add the segment
      fseed._segments.emplace_back(zpiece,tz); //tmin,tmax
    }
    return fseed;
  }

  //TODO - the followinf function requires work:
  // 1) KinematicLine - constructor doesnt work
  // 2) paramVal-->paramVar
  // 3) set parameters
  KTRAJ KinematicLineFit::makeSeedTraj(CosmicTrackSeed const& hseed) const { //TODO - switch to CosmicTrackSeed
    // first, find the time range of the hits
    auto const& hhits = hseed.hits();
    float tmin = std::numeric_limits<float>::max();
    float tmax = std::numeric_limits<float>::min();
    for( auto const& hhit: hhits) {
      tmin = std::min(tmin,hhit.correctedTime());
      tmax = std::max(tmax,hhit.correctedTime());
    }
    TimeRange trange(tmin-tbuff_,tmax+tbuff_);
    //auto const& scosmic = hseed.track(); //a CosmicTrack

    VEC3 bnom(0.0,0.0,0.0);
    // create a PKTRAJ from the CosmicTrack fit result, to seed the KinKal fit.  First, translate the parameters
    std::tuple <double, double, double, double, double> info = scosmic.GetTrackPOCAInfo();//d0,phi0,z0,cost,t0
    DVEC pars;
    pars[KTRAJ::d0_] = get<0>(info);
    pars[KTRAJ::z0_] = get<2>(info);
    pars[KTRAJ::cost_] = get<3>(info);
    pars[KTRAJ::phi0_] = get<1>(info);
    pars[KTRAJ::mom_] = 1.0; //TODO
    pars[KTRAJ::t0_] = get<4>(info); //TODO

    // create the initial trajectory
    Parameters kkpars(pars,seedcov_);
    if(print_ > 1){
      std::cout << "Seed Cosmic parameters " << scosmic << std::endl;
      std::cout << "Seed Traj parameters  " << kkpars << std::endl;
    }
    //  construct the seed trajectory
    return KTRAJ(kkpars, bnom, 0.0 );//TODO
  }

  //Following functions are the same as for LoopHelix:
  void KinematicLineFit::makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,KTRAJ const& seedtraj,
      ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
      MEASCOL& thits,EXINGCOL& exings) const {
        // wrap the seed traj in a Piecewise traj: needed to satisfy PTOCA interface
        PKTRAJ pseedtraj(seedtraj);
        // loop over the individual straw hits
        for(auto strawidx : strawHitIdxs) {
          const ComboHit& strawhit(chcol.at(strawidx));
          if(strawhit.mask().level() != StrawIdMask::uniquestraw)
          throw cet::exception("RECO")<<"mu2e::KinematicLineFit: ComboHit error"<< endl;
          const Straw& straw = tracker.getStraw(strawhit.strawId());
          // find the propagation speed for signals along this straw
          double sprop = 2*strawresponse.halfPropV(strawhit.strawId(),strawhit.energyDep());
          // construct a kinematic line trajectory from this straw. the measurement point is the signal end
          auto p0 = straw.wireEnd(strawhit.driftEnd());
          auto p1 = straw.wireEnd(StrawEnd(strawhit.driftEnd().otherEnd()));
          auto propdir = (p0 - p1).unit(); // The signal propagates from the far end to the near
          // clumsy conversion: make this more elegant TODO
          VEC3 vp0(p0.x(),p0.y(),p0.z());
          VEC3 vp1(p1.x(),p1.y(),p1.z());
          Line wline(vp0,vp1,strawhit.time(),sprop);
          // compute 'hint' to TPOCA.  correct the hit time using the time division
          double psign = propdir.dot(straw.wireDirection());  // wire distance is WRT straw center, in the nominal wire direction
          double htime = wline.t0() - (straw.halfLength()-psign*strawhit.wireDist())/wline.speed();
          CAHint hint(seedtraj.ztime(vp0.Z()),htime);
          // compute a preliminary PTCA between the seed trajectory and this straw.
          PTCA ptca(pseedtraj, wline, hint, tpocaprec_);
          // create the material crossing
          if(addmat_){
            exings.push_back(std::make_shared<STRAWXING>(ptca,*smat_));
            if(print_ > 2)exings.back()->print(std::cout,1);
          }
          // create the hit
          thits.push_back(std::make_shared<KKSTRAWHIT>(*kkbf_, ptca, whstate_, strawhit, straw, strawidx, strawresponse));
          if(print_ > 2)thits.back()->print(std::cout,2);
        }
  }

  double KinematicLineFit::zTime(PKTRAJ const& ptraj, double zpos) const {
    auto bpos = ptraj.position3(ptraj.range().begin());
    auto epos = ptraj.position3(ptraj.range().end());
    // assume linear transit to get an initial estimate
    double tz = ptraj.range().begin() + ptraj.range().range()*(zpos-bpos.Z())/(epos.Z()-bpos.Z());
    size_t zindex = ptraj.nearestIndex(tz);
    auto const& traj = ptraj.piece(zindex);
    bpos = traj.position3(traj.range().begin());
    epos = traj.position3(traj.range().end());
    tz = traj.range().begin() + traj.range().range()*(zpos-bpos.Z())/(epos.Z()-bpos.Z());
    return tz;
  }


}
// mu2e

DEFINE_ART_MODULE(mu2e::KinematicLineFit);
