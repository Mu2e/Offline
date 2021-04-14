#ifndef Mu2eKinKal_KKFit_hh
#define Mu2eKinKal_KKFit_hh
//
// helper class for constructing KinKal Fits using Mu2e data
//
#include "Mu2eKinKal/inc/KKStrawHit.hh"
//#include "Mu2eKinKal/inc/KKPanelHit.hh"
#include "Mu2eKinKal/inc/KKCaloHit.hh"
#include "Mu2eKinKal/inc/KKFitSettings.hh"
#include "Mu2eKinKal/inc/KKBField.hh"
// mu2e Includes
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
// KinKal includes
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/Status.hh"
#include "KinKal/Fit/Track.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/StrawXing.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
namespace mu2e {
  using KinKal::WireHitState;
  using KinKal::Line;
  using KinKal::StrawMaterial;
  using KinKal::TimeRange;
  using KinKal::Config;
  using KinKal::Status;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;
  using KKFitSettings = Mu2eKinKal::FitSettings;

  template <class KTRAJ> class KKFit {
    public:
      // fit configuration
      using KKTRK = KinKal::Track<KTRAJ>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using TCA = KinKal::ClosestApproach<KTRAJ,Line>;
      using KKHIT = KinKal::HitConstraint<KTRAJ>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKCALOHIT = KKCaloHit<KTRAJ>;
      //      using KKPANELHIT = KKPanelHit<KTRAJ>;
      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using STRAWXING = KinKal::StrawXing<KTRAJ>;
      using HPtr = art::Ptr<HelixSeed>;
      // construct from a fit configuration objects
      explicit KKFit(KKFitSettings const& fitconfig);
      // helper functions used to create components of the fit
      void makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,KKBField const& kkbf, StrawMaterial const& smat,
	  PKTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
	  MEASCOL& hits, EXINGCOL& exings) const;
      void makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PKTRAJ const& pktraj, MEASCOL& hits) const;
      KalSeed createSeed(KKTRK const& kktrk, HPtr const& hptr, std::vector<float> const& zsave, bool savefull) const;
      TimeRange range(MEASCOL const& hits, EXINGCOL const& xings) const; // time range from a set of hits and element Xings
      double zTime(PKTRAJ const& trak, double zpos) const; // find the time the trajectory crosses the plane perp to z at the given z position
      bool useCalo() const { return usecalo_; }
      WireHitState::Dimension nullDimension() const { return nulldim_; }
      Config const& config() const { return config_; }
      Config& config() { return config_; }
      PDGCode::type fitParticle() const { return tpart_;}
      TrkFitDirection fitDirection() const { return tdir_;}
    private:
      PDGCode::type tpart_;
      TrkFitDirection tdir_;
      WireHitState::Dimension nulldim_;
      float nullvscale_;
      bool addmat_, usecalo_;
      double caloDt_; // calo time offset; should come from proditions FIXME!
      double caloPosRes_; // calo cluster transverse position resolution; should come from proditions or CaloCluster FIXME!
      double caloPropSpeed_; // effective light propagation speed in a crystal (including reflections).  Should come from prodtions FIXME
      Config config_; // configuration object
  };

  template <class KTRAJ> KKFit<KTRAJ>::KKFit(KKFitSettings const& fitconfig) :
    tpart_(static_cast<PDGCode::type>(fitconfig.fitParticle())),
    tdir_(static_cast<TrkFitDirection::FitDirection>(fitconfig.fitDirection())),
    nulldim_(static_cast<WireHitState::Dimension>(fitconfig.nullHitDimension())),
    nullvscale_(fitconfig.nullHitVarianceScale()),
    addmat_(fitconfig.addMaterial()),
    usecalo_(fitconfig.useCaloCluster()),
    caloDt_(fitconfig.caloDt()),
    caloPosRes_(fitconfig.caloPosRes()),
    caloPropSpeed_(fitconfig.caloPropSpeed())
  {
    // fill configuration
    config_.maxniter_ = fitconfig.maxniter();
    config_.tprec_ = fitconfig.tpocaPrec();
    config_.dwt_ = fitconfig.dwt();
    config_.pdchi2_ = fitconfig.dparams();
    config_.tbuff_ = fitconfig.tBuffer();
    config_.tol_ = fitconfig.btol();
    config_.minndof_ = fitconfig.minndof();
    config_.bfcorr_ = static_cast<Config::BFCorr>(fitconfig.bfieldCorr());
    config_.plevel_ = static_cast<Config::printLevel>(fitconfig.printLevel());

    // set the schedule for the meta-iterations
    unsigned nmiter(0);
    for(auto const& misetting : fitconfig.mconfig()) {
      MetaIterConfig mconfig;
      mconfig.temp_ = std::get<0>(misetting);
      mconfig.convdchisq_ = std::get<1>(misetting);
      mconfig.divdchisq_ = std::get<2>(misetting);
      mconfig.miter_ = nmiter++;
      config_.schedule_.push_back(mconfig);
    }
  }

  template <class KTRAJ> void KKFit<KTRAJ>::makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,KKBField const& kkbf, StrawMaterial const& smat,
      PKTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
      MEASCOL& hits, EXINGCOL& exings) const {
    // initialize hits as null (no drift).  Drift is turned on when updating
    auto const& sprop = tracker.strawProperties();
    double rstraw = sprop.strawInnerRadius();
    double nulldt = 0.5*rstraw/strawresponse.driftConstantSpeed(); // approximate shift in time due to ignoring drift
    double nullvar = nullvscale_*rstraw*rstraw/3.0; // scaled square RMS (distance is between 0 and r)
    WireHitState whstate(WireHitState::null, nulldim_, nullvar ,nulldt); // initial wire hit state

    // loop over the individual straw hits
    for(auto strawidx : strawHitIdxs) {
      const ComboHit& strawhit(chcol.at(strawidx));
      if(strawhit.mask().level() != StrawIdMask::uniquestraw)
	throw cet::exception("RECO")<<"mu2e::KKFit: ComboHit error"<< endl;
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
      CAHint hint(ptraj.front().ztime(vp0.Z()),htime);
      // compute a preliminary PTCA between the seed trajectory and this straw.
      PTCA ptca(ptraj, wline, hint, config_.tprec_ );
      // create the material crossing
      if(addmat_){
	exings.push_back(std::make_shared<STRAWXING>(ptca,smat));
      }
      // create the hit
      hits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, ptca, whstate, strawhit, straw, strawidx, strawresponse));
    }
  }

  template <class KTRAJ> void KKFit<KTRAJ>::makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PKTRAJ const& ptraj, MEASCOL& hits) const {
    // move cluster COG into the tracker frame.  COG is at the front face of the disk
    CLHEP::Hep3Vector cog = calo.geomUtil().mu2eToTracker(calo.geomUtil().diskFFToMu2e( cluster->diskID(), cluster->cog3Vector()));
    // project this along the crystal axis to the SIPM, which is at the back.  This is the point the time measurement corresponds to
    VEC3 FFCOG(cog);
    double lcrystal = calo.caloInfo().getDouble("crystalZLength");
    VEC3 crystalF2B = VEC3(0.0,0.0,lcrystal); // this should come directly from the calogeometry, TODO
    VEC3 SIPMCOG = FFCOG + crystalF2B;
    // create the Line trajectory from this information: signal goes towards the sipm
    Line caxis(SIPMCOG,FFCOG,cluster->time()+caloDt_,caloPropSpeed_); 
    CAHint hint( caxis.t0(), caxis.t0());
    // compute a preliminary PTCA between the seed trajectory and this straw.
    PTCA ptca(ptraj, caxis, hint, config_.tprec_ );
    // create the hit
    double tvar = cluster->timeErr()*cluster->timeErr();
    double wvar = caloPosRes_*caloPosRes_;
    hits.push_back(std::make_shared<KKCALOHIT>(cluster,ptca,tvar,wvar));
  }

  template <class KTRAJ> double KKFit<KTRAJ>::zTime(PKTRAJ const& ptraj, double zpos) const {
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

  template <class KTRAJ> TimeRange KKFit<KTRAJ>::range(MEASCOL const& hits, EXINGCOL const& xings) const{
    double tmin = std::numeric_limits<float>::max();
    double tmax = std::numeric_limits<float>::min();
    for( auto const& hit : hits) {
      tmin = std::min(tmin,hit->time());
      tmax = std::max(tmax,hit->time());
    }
    for( auto const& xing : xings) {
      tmin = std::min(tmin,xing->crossingTime());
      tmax = std::max(tmax,xing->crossingTime());
    }
    return TimeRange(tmin-config_.tbuff_,tmax+config_.tbuff_);
  }

  template <class KTRAJ> KalSeed KKFit<KTRAJ>::createSeed(KKTRK const& kktrk,HPtr const& hptr, std::vector<float> const& zsave, bool savefull) const {
    TrkFitFlag fflag(hptr->status());
    if(kktrk.fitStatus().usable()) fflag.merge(TrkFitFlag::kalmanOK);
    if(kktrk.fitStatus().status_ == Status::converged) fflag.merge(TrkFitFlag::kalmanConverged);
    if(addmat_)fflag.merge(TrkFitFlag::MatCorr);
    if(config().bfcorr_ != 0)fflag.merge(TrkFitFlag::BFCorr);
    // explicit T0 is needed for backwards-compatibility; sample from the appropriate trajectory piece
    auto const& fittraj = kktrk.fitTraj();
    double tz0 = zTime(fittraj,0.0);
    auto const& t0piece = fittraj.nearestPiece(tz0);
    HitT0 t0(t0piece.paramVal(KTRAJ::t0_), sqrt(t0piece.paramVar(KTRAJ::t0_))); 
    // create the shell for the output.  Note the (obsolete) flight length is given as t0
    KalSeed fseed(tpart_,tdir_,fflag,t0.t0());
    fseed._helix = hptr;
    auto const& fstatus = kktrk.fitStatus();
    fseed._chisq = fstatus.chisq_.chisq();
    fseed._fitcon = fstatus.chisq_.probability();
    fseed._nseg = fittraj.pieces().size();
    // loop over individual effects
    for(auto const& eff: kktrk.effects()) {
      const KKHIT* kkhit = dynamic_cast<const KKHIT*>(eff.get());
      if(kkhit != 0){
	const KKSTRAWHIT* strawhit = dynamic_cast<const KKSTRAWHIT*>(kkhit->hit().get());
	const KKCALOHIT* calohit = dynamic_cast<const KKCALOHIT*> (kkhit->hit().get());
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
	} else if(calohit != 0) {
	  auto const& ca = calohit->closestApproach();
	  StrawHitFlag hflag;
	  if(calohit->active()){
	    hflag.merge(StrawHitFlag::active);
	    hflag.merge(StrawHitFlag::doca);
	  }
	  fseed._chit = TrkCaloHitSeed(HitT0(ca.particleToca(),sqrt(ca.tocaVar())),
	      static_cast<float>(ca.particleToca()), static_cast<float>(ca.sensorToca()),
	      static_cast<float>(ca.doca()),static_cast<float>(sqrt(ca.docaVar())), static_cast<float>(ca.sensorToca()),
	      static_cast<float>(sqrt(ca.tocaVar())),hflag);
	  fseed._chit._cluster = calohit->caloCluster();
	  // save the segment at the CH
	  auto const& zpiece = fittraj.nearestPiece(ca.particleToca());
	  fseed._segments.emplace_back(zpiece,ca.particleToca());
	}
      }
      //      const KKMAT* kkmat = dynamic_cast<const KKMAT*>(eff.get()); TODO!!
      //TrkUtilities::fillStraws(ktrk,fseed._straws); TODO!!
    }
    // sample the fit at the requested z positions.
    if(savefull){
      // loop over all pieces of the fit trajectory
      for (auto const& traj : fittraj.pieces() ) {
	fseed._segments.emplace_back(traj,traj.range().mid());
      }
    } else {
      for(auto zpos : zsave ) {
	// compute the time the trajectory crosses this plane
	double tz = zTime(fittraj,zpos);
	// find the explicit trajectory piece at this time
	auto const& zpiece = fittraj.nearestPiece(tz);
	// construct and add the segment
	fseed._segments.emplace_back(zpiece,tz);
      }
    }
    return fseed;
  }

}
#endif
