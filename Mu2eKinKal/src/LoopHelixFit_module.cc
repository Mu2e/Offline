//
// KinKal fit module using the LoopHelix parameterset
//
// Original author D. Brown (LBNL) 11/18/20
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
#include "KinKal/Detector/StrawMat.hh"
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/Trajectory/Line.hh"
// Mu2eKinKal
#include "Mu2eKinKal/inc/KKFileFinder.hh"
#include "Mu2eKinKal/inc/KKStrawHit.hh"
//#include "Mu2eKinKal/inc/KKPanelHit.hh"
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
//  using TSH = KKStrawHit<KTRAJ>;
//  using TPH = KKPanelHit<KTRAJ>;
  using KKTRK = KinKal::Track<KTRAJ>;
  using STRAWXING = KinKal::StrawXing<KTRAJ>;
  using KinKal::Line;
  using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
  using MatEnv::MatDBInfo;
  using KinKal::DVEC;
  using KinKal::Parameters;
  using KinKal::VEC3;
  using KinKal::TimeRange;
  using KinKal::Config;
  using KinKal::MetaIterConfig;
  using KinKal::CAHint;
  using KinKal::StrawMat;
  using MEAS = KinKal::Hit<KTRAJ>;
  using MEASPTR = std::shared_ptr<MEAS>;
  using MEASCOL = std::vector<MEASPTR>;
  using KKSTRAWHIT = KKStrawHit<KTRAJ>;

  class LoopHelixFit : public art::EDProducer {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    struct ModuleSettings {
      fhicl::Sequence<art::InputTag> seedCollections         {Name("HelixSeedCollections"),     Comment("Helix seed fit collections to be processed ") };
      fhicl::Atom<art::InputTag>     comboHitCollection     {Name("ComboHitCollection"),     Comment("Single Straw ComboHit collection ") };
      fhicl::Atom<art::InputTag>     strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection ") };
      fhicl::Atom<int> fitParticle {  Name("fitparticle"), Comment("Particle type to fit: e-, e+, mu-, ..."), PDGCode::e_minus};
      fhicl::Atom<int> fitDirection { Name("fitdirection"), Comment("Particle direction to fit, either upstream or downstream"), TrkFitDirection::downstream };
      fhicl::Atom<bool> refine { Name("Refine"), Comment("Refine and redo the final fit") };
      fhicl::Atom<bool> useCalo { Name("UseCaloCluster"), Comment("Add CaloCluster info to final fit or not") };
      fhicl::Atom<int> diagLevel { Name("DiagLevel"), Comment("Diagnostic Level"), 0 };
      fhicl::Atom<int> debugLevel { Name("DebugLevel"), Comment("Debug Level"), 0 };
      fhicl::Sequence<std::string> helixFlags { Name("HelixFlags"), Comment("Flags required to be present to convert a seed to a KinKal track") };
      fhicl::Atom<bool> saveAll { Name("SaveAllFits"), Comment("Save all fits, whether they suceed or not"),false };
      fhicl::Sequence<std::string> addHitFlags { Name("AddHitFlags"), Comment("Flags required to be present to add a hit") };
      fhicl::Sequence<std::string> rejectHitFlags { Name("RejectHitFlags"), Comment("Flags required not to be present to add a hit") };
      fhicl::Atom<float> maxAddDOCA { Name("MaxAddDOCA"), Comment("Max DOCA to add a hit") };
      fhicl::Atom<float> maxAddDt { Name("MaxAddDt"), Comment("Max Detla time to add a hit") };
      fhicl::Atom<float> maxAddChi { Name("MaxAddChi"), Comment("Max Chi to add a hit") };
      fhicl::Atom<float> maxAddDeltaU { Name("MaxAddDeltaU"), Comment("Max Delta-U to add a hit") };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for initiaion time range") };
      fhicl::Atom<float> tPOCAPrec { Name("TPOCAPrecision"), Comment("Precision for initial TPOCA (ns)"), 0.1 };
      fhicl::Sequence<float> zsave { Name("ZSavePositions"), Comment("Z positions to sample and save the fit")};
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
      fhicl::Atom<int> maxniter { Name("MaxNIter"), Comment("Maximum number of algebraic iteration steps in each fit meta-iteration"), 10 };
      fhicl::Atom<float> dwt { Name("Deweight"), Comment("Deweighting factor when initializing the track end parameters"), 1.0e6 };
      fhicl::Atom<float> tBuffer { Name("TimeBuffer"), Comment("Time buffer for final fit") };
      fhicl::Atom<float> btol { Name("BCorrTolerance"), Comment("Tolerance on BField correction accuracy (mm)"), 0.01 };
      fhicl::Atom<int> minndof { Name("MinNDOF"), Comment("Minimum number of Degrees of Freedom to conitnue fitting"), 5  };
      fhicl::Atom<bool> addmat { Name("AddMaterial"), Comment("Add material effecst to the fit"), true };
      fhicl::Atom<int> bfieldCorr { Name("BFieldCorr"), Comment("BField correction algorith") };
      fhicl::Atom<float> nullvar { Name("NullHitVariance"), Comment("Spatial variance for null wire hits (mm)"), 1.0 };
      fhicl::Atom<float> nulltoffset { Name("NullHitTimeOffset"), Comment("Time Offset for null wire hits (ns)"), 20.0 };
      fhicl::Atom<int> printLevel { Name("PrintLevel"), Comment("Print Level"),0};
    };

    using MetaSettings = fhicl::Sequence<fhicl::Tuple<bool,bool,bool,float,float,float,float,float>>;
    struct ModuleConfig {
      fhicl::Table<ModuleSettings> modsettings { Name("ModuleSettings") };
      fhicl::Table<FitSettings> fitsettings { Name("FitSettings") };
      fhicl::Table<MaterialSettings> matsettings { Name("MaterialSettings") };
      MetaSettings mconfig { Name("MetaIterationSettings"), Comment("MetaIteration sequence configuration parameters, format: \n"
      " 'UpdateMaterials (t/f)', 'UpdateBFieldCorrection (t/f)', 'UpdateHits (t/f)', 'Temperature (dimensionless)', \n"
      "'chisq/dof for convergence', 'chisq/dof for divergence', 'chisq/dof for oscillation', 'TPOCA convergence tolerance (ns)'") };
    };
    using ModuleParams = art::EDProducer::Table<ModuleConfig>;

    public:
    explicit LoopHelixFit(const ModuleParams& config);
    virtual ~LoopHelixFit();
    void beginRun(art::Run& aRun) override;
    void produce(art::Event& event) override;
    private:
    std::vector<art::ProductToken<HelixSeedCollection>> seedCols_;
    art::ProductToken<ComboHitCollection> chcol_T_;
    art::ProductToken<StrawHitFlagCollection> shfcol_T_;
    TrkFitFlag goodhelix_;
    bool saveall_;
    TrkFitDirection tdir_;
    PDGCode::type tpart_;
    ProditionsHandle<StrawResponse> strawResponse_h_;
    ProditionsHandle<Tracker> alignedTracker_h_;
    float maxDoca_, maxDt_, maxChi_, maxDU_, tbuff_, tpocaprec_;
    KKFileFinder filefinder_;
    std::string wallmatname_, gasmatname_, wirematname_;
    std::unique_ptr<StrawMat> smat_; // straw material
    Config config_;
    KinKal::WireHitState whstate_; // initial state of wire hits
  };

  LoopHelixFit::LoopHelixFit(const ModuleParams& config) : art::EDProducer{config}, 
    chcol_T_(consumes<ComboHitCollection>(config().modsettings().comboHitCollection())),
    shfcol_T_(mayConsume<StrawHitFlagCollection>(config().modsettings().strawHitFlagCollection())),
    goodhelix_(config().modsettings().helixFlags()),
    saveall_(config().modsettings().saveAll()),
    tdir_(static_cast<TrkFitDirection::FitDirection>(config().modsettings().fitDirection())), tpart_(static_cast<PDGCode::type>(config().modsettings().fitParticle())),
    maxDoca_(config().modsettings().maxAddDOCA()),
    maxDt_(config().modsettings().maxAddDt()),
    maxChi_(config().modsettings().maxAddChi()),
    maxDU_(config().modsettings().maxAddDeltaU()),
    tbuff_(config().modsettings().tBuffer()),
    tpocaprec_(config().modsettings().tPOCAPrec()),
    filefinder_(config().matsettings().elements(),config().matsettings().isotopes(),config().matsettings().materials()),
    wallmatname_(config().matsettings().strawWallMaterialName()),
    gasmatname_(config().matsettings().strawGasMaterialName()),
    wirematname_(config().matsettings().strawWireMaterialName()),
    whstate_(KinKal::WireHitState::null, KinKal::WireHitState::both,config().fitsettings().nullvar(), config().fitsettings().nulltoffset())
  {
    // collection handling
    for(const auto& seedtag : config().modsettings().seedCollections()) { seedCols_.emplace_back(consumes<HelixSeedCollection>(seedtag)); }
    produces<KKLoopHelixCollection>();
    produces<KalSeedCollection>();
    // construct the fit configuration object.  This controls all the global and iteration-specific aspects of the fit
    config_.maxniter_ = config().fitsettings().maxniter();
    config_.dwt_ = config().fitsettings().dwt();
    config_.tbuff_ = config().fitsettings().tBuffer();
    config_.tol_ = config().fitsettings().btol();
    config_.minndof_ = config().fitsettings().minndof();
    config_.addmat_ = config().fitsettings().addmat();
    config_.bfcorr_ = static_cast<Config::BFCorr>(config().fitsettings().bfieldCorr());
    config_.plevel_ = static_cast<Config::printLevel>(config().fitsettings().printLevel());
    // Now set the schedule for the meta-iterations
    unsigned nmiter(0);
    for(auto const& misetting : config().mconfig()) {
      MetaIterConfig mconfig;
      mconfig.updatemat_ = std::get<0>(misetting);
      mconfig.updatebfcorr_ = std::get<1>(misetting);
      mconfig.updatehits_ = std::get<2>(misetting);
      mconfig.temp_ = std::get<3>(misetting);
      mconfig.convdchisq_ = std::get<4>(misetting);
      mconfig.divdchisq_ = std::get<5>(misetting);
      mconfig.oscdchisq_ = std::get<6>(misetting);
      mconfig.tprec_ = std::get<7>(misetting);
      mconfig.miter_ = nmiter++;
      config_.schedule_.push_back(mconfig);
    }
  }

  LoopHelixFit::~LoopHelixFit(){}
  //-----------------------------------------------------------------------------
  void LoopHelixFit::beginRun(art::Run& run) {
    // initialize material access
    GeomHandle<Tracker> tracker_h;
    auto const& sprop = tracker_h->strawProperties();
    MatDBInfo matdbinfo(filefinder_);
    smat_ = std::make_unique<StrawMat>(
	sprop._strawOuterRadius, sprop._strawWallThickness, sprop._wireRadius,
	matdbinfo.findDetMaterial(wallmatname_),
	matdbinfo.findDetMaterial(gasmatname_),
	matdbinfo.findDetMaterial(wirematname_));
  }

  void LoopHelixFit::produce(art::Event& event ) {
    // find current proditions
    auto const& ptable = GlobalConstantsHandle<ParticleDataTable>();
    auto const& strawresponse = strawResponse_h_.getPtr(event.id());
    auto const& tracker = alignedTracker_h_.getPtr(event.id()).get();
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    KKBField kkbf(*bfmgr,*det); // should set this at beginrun TODO
    // find input hits
    auto ch_H = event.getValidHandle<ComboHitCollection>(chcol_T_);
    auto const& chcol = *ch_H;;

    // create output
    unique_ptr<KKLoopHelixCollection> kktrkcol(new KKLoopHelixCollection );
    unique_ptr<KalSeedCollection> kkseedcol(new KalSeedCollection );
    // find the seed collections
    for (auto const& seedtag : seedCols_) {
      auto const& seedcol_h = event.getValidHandle<HelixSeedCollection>(seedtag);
      // loop over the seeds
      for( auto const& seed : *seedcol_h ) {
	if(seed.status().hasAllProperties(goodhelix_)){
	  // create a PKTRAJ from the helix fit result, to seed the KinKal fit.  First, translate the parameters
	  // Note the flips in case of backwards propagation
	  auto const& shelix = seed.helix();
	  DVEC pars;
	  pars[KTRAJ::rad_] = shelix.radius();
	  pars[KTRAJ::lam_] = shelix.lambda();
	  pars[KTRAJ::cx_] = shelix.centerx();
	  pars[KTRAJ::cy_] = shelix.centery();
	  pars[KTRAJ::phi0_] = shelix.fz0();
	  pars[KTRAJ::t0_] = seed.t0().t0();
	  if(tdir_ == TrkFitDirection::upstream){
	    pars[KTRAJ::rad_] *= -1.0;
	    pars[KTRAJ::lam_] *= -1.0;
	  }
	  Parameters kkpars(pars);
	  // compute the magnetic field at the center.  We only want the z compontent, as the helix fit assumes B points along Z
	  float zcent = 0.5*(seed.hits().front().pos().Z()+seed.hits().back().pos().Z());
	  VEC3 center(shelix.centerx(), shelix.centery(),zcent);
	  auto bcent = kkbf.fieldVect(center);
	  VEC3 bnom(0.0,0.0,bcent.Z());
	  // now compute kinematics
	  double mass = ptable->particle(tpart_).ref().mass().value(); 
	  int charge = static_cast<int>(ptable->particle(tpart_).ref().charge());
	  // find the time range, using the hit time
	  //	float tmin = std::numeric_limits<float>::max();
	  //	float tmax = std::numeric_limits<float>::min();
	  TimeRange trange(seed.hits().front().correctedTime() - tbuff_, seed.hits().back().correctedTime() + tbuff_);
	  //  construct the trajectory
	  KTRAJ ktraj(kkpars, trange, mass, charge, bnom);
	  PKTRAJ seedtraj(ktraj);
	  // construct individual KKHit objects from each Helix hit
	  // first, we need to unwind the straw combination
	  std::vector<StrawHitIndex> strawHitIdxs;
	  for(size_t ihit = 0; ihit < seed.hits().size(); ++ihit ){
// dig down to the straw-level hits
	    seed.hits().fillStrawHitIndices(event,ihit,strawHitIdxs);
	  }
	  MEASCOL thits; //polymorphic container of hits

	  for(auto strawidx : strawHitIdxs) {
	    const ComboHit& strawhit(chcol.at(strawidx));
	    if(strawhit.mask().level() != StrawIdMask::uniquestraw)
	      throw cet::exception("RECO")<<"mu2e::LoopHelixFit: ComboHit error"<< endl;
	    const Straw& straw = tracker->getStraw(strawhit.strawId());
	    // find the propagation speed for signals along this straw
	    double sprop = 2*strawresponse->halfPropV(strawhit.strawId(),strawhit.energyDep());
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
	    CAHint hint(ktraj.ztime(vp0.Z()),htime);
	    // compute a preliminary PTCA between the seed trajectory and this straw.
	    PTCA ptca(seedtraj, wline, hint, tpocaprec_);
	    std::cout << ptca.doca() << std::endl;
	    // create the material crossing from this PTCA
	    auto sxing = std::make_shared<STRAWXING>(ptca,*smat_);
	    // now create the hit; initial state is 
	    thits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, ptca, sxing, whstate_, strawhit, straw, *strawresponse));
//		  sdrift_, sigt_*sigt_, rstraw_));

	  }
	  // search for straws near the trajectory, and add them as passive material if they aren't already hits


	  // build the KinKal track.  Fit is performed on construction
	  //	  kktrk = new KKTRACK(config_, 
	  //	      art::Ptr<CaloCluster> ccPtr;
	  //	      if (kseed.caloCluster()){
	  //	      ccPtr = kseed.caloCluster(); // remember the Ptr for creating the TrkCaloHitSeed and KalSeed Ptr
	  //	      }
	  //
	  //	      // save all the fits
	  //	      kktrkcol->push_back(kktrk);
	  //
	  //	      // convert successful fits into 'seeds' for persistence
	  //	      TrkFitFlag fflag(kseed.status());
	  //	      fflag.merge(TrkFitFlag::KFF);
	  //	      KalSeed fseed(krep->particleType(),_fdir,krep->t0(),krep->flt0(),fflag);
	  //
	  //
	  //	      // save KalSeed for this track
	  //	      kscol->push_back(fseed)
	  //	      if (diag_ > 0) _hmanager->fillHistograms(&_data);
	}
      }
    }

    // put the output products into the event
    //  event.put(move(krcol));
  }

  //
  //void LoopHelixFit::createKalSeed(KKTRACK const& kktrk, KalSeed& kseed) {
  //  if(krep->fitStatus().success()) fflag.merge(TrkFitFlag::kalmanOK);
  //  if(krep->fitStatus().success()==1) fflag.merge(TrkFitFlag::kalmanConverged);
  //
  //
  //  // reference the seed fit in this fit
  //  auto ksH = event.getValidHandle<KalSeedCollection>(_ksToken);
  //  fseed._kal = art::Ptr<KalSeed>(ksH,ikseed);
  //  // redundant but possibly useful
  //  fseed._helix = kseed.helix();
  //  // fill with new information
  //  fseed._t0 = krep->t0();
  //  fseed._flt0 = krep->flt0();
  //  // global fit information
  //	    fseed._chisq = krep->chisq();
  //	    // compute the fit consistency.  Note our fit has effectively 6 parameters as t0 is allowed to float and its error is propagated to the chisquared
  //	    fseed._fitcon =  TrkUtilities::chisqConsistency(krep);
  //	    fseed._nbend = TrkUtilities::countBends(krep);
  //	    TrkUtilities::fillStrawHitSeeds(krep,*_chcol,fseed._hits);
  //	    TrkUtilities::fillStraws(krep,fseed._straws);
  //	    // sample the fit at the requested z positions.  Need options here to define a set of
  //	    // standard points, or to sample each unique segment on the fit FIXME!
  //	    for(auto zpos : _zsave) {
  //	      // compute the flightlength for this z
  //	      double fltlen = krep->pieceTraj().zFlight(zpos);
  //	      // sample the momentum at this flight.  This belongs in a separate utility FIXME
  //	      BbrVectorErr momerr = krep->momentumErr(fltlen);
  //	      // sample the helix
  //	      double locflt(0.0);
  //	      const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(fltlen,locflt));
  //	      // fill the segment
  //	      KalSegment kseg;
  //	      TrkUtilities::fillSegment(*htraj,momerr,locflt-fltlen,kseg);
  //	      fseed._segments.push_back(kseg);
  //	    }
  //	    // see if there's a TrkCaloHit
  //	    const TrkCaloHit* tch = TrkUtilities::findTrkCaloHit(krep);
  //	    if(tch != 0){
  //	      TrkUtilities::fillCaloHitSeed(tch,fseed._chit);
  //	      // set the Ptr using the helix: this could be more direct FIXME!
  //	      fseed._chit._cluster = ccPtr;
  //	      // create a helix segment at the TrkCaloHit
  //	      KalSegment kseg;
  //	      // sample the momentum at this flight.  This belongs in a separate utility FIXME
  //	      BbrVectorErr momerr = krep->momentumErr(tch->fltLen());
  //	      double locflt(0.0);
  //	      const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(tch->fltLen(),locflt));
  //	      TrkUtilities::fillSegment(*htraj,momerr,locflt-tch->fltLen(),kseg);
  //	      fseed._segments.push_back(kseg);
  //	    }
  //
  //  }
  //
  //  void LoopHelixFit::findMissingHits(KalFitData&kalData) {
  //    KinKalTrack* krep = kalData.krep;
//
//    //clear the array
//    kalData.missingHits.clear();
//    const Tracker& tracker = *_data.tracker;
//
//    //  Trajectory info
//    Hep3Vector tdir;
//    HepPoint tpos;
//    krep->pieceTraj().getInfo(krep->flt0(),tpos,tdir);
//    unsigned nstrs = _chcol->size();
//    TrkStrawHitVector tshv;
//    convert(krep->hitVector(),tshv);
//    for(unsigned istr=0; istr<nstrs;++istr){
//      if(_shfcol->at(istr).hasAllProperties(_addsel)&& !_shfcol->at(istr).hasAnyProperty(_addbkg)){
//	ComboHit const& sh = _chcol->at(istr);
//	if (sh.flag().hasAnyProperty(StrawHitFlag::dead)) {
//	  continue;
//	}
//	if(fabs(_chcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
//	  // make sure we haven't already used this hit
//	  vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
//	  if(ifnd == tshv.end()){
//	    // good in-time hit.  Compute DOCA of the wire to the trajectory
//	    Straw const& straw = tracker.getStraw(sh.strawId());
//	    CLHEP::Hep3Vector hpos = straw.getMidPoint();
//	    CLHEP::Hep3Vector hdir = straw.getDirection();
//	    // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
//	    HepPoint spt(hpos.x(),hpos.y(),hpos.z());
//	    TrkLineTraj htraj(spt,hdir,-straw.halfLength(),straw.halfLength());
//	    // estimate flightlength along track.  This assumes a constant BField!!!
//	    double fltlen = (hpos.z()-tpos.z())/tdir.z();
//	    // estimate hit length
//	    HepPoint tp = krep->pieceTraj().position(fltlen);
//	    Hep3Vector tpos(tp.x(),tp.y(),tp.z()); // ugly conversion FIXME!
//	    double hitlen = hdir.dot(tpos - hpos);
//	    TrkPoca hitpoca(krep->pieceTraj(),fltlen,htraj,hitlen);
//
//	    // flag hits with small residuals
//	    if(fabs(hitpoca.doca()) < _maxadddoca){
//	      MissingHit_t m;
//	      m.index = istr;
//	      m.doca  = hitpoca.doca();
//	      // m.dr = ??;
//	      kalData.missingHits.push_back(m);
//	    }
//	  }
//	}
//      }
//    }
//  }
//

}
// mu2e

DEFINE_ART_MODULE(mu2e::LoopHelixFit);
