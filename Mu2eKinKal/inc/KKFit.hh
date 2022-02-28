#ifndef Mu2eKinKal_KKFit_hh
#define Mu2eKinKal_KKFit_hh
//
// helper class for constructing KinKal Fits using Mu2e data
//
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
//#include "Mu2eKinKal/inc/KKPanelHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/KKCaloHit.hh"
#include "Offline/Mu2eKinKal/inc/KKFitUtilities.hh"
#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKBField.hh"
// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"
// mu2e Includes
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
// KinKal includes
#include "KinKal/Fit/Status.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Detector/StrawMaterial.hh"
// Other
#include "cetlib_except/exception.h"
#include <memory>
#include <cmath>
#include <algorithm>
namespace mu2e {
  using KinKal::WireHitState;
  using KinKal::Line;
  using KinKal::StrawMaterial;
  using KinKal::TimeRange;
  using KinKal::Status;
  using StrawHitIndexCollection = std::vector<StrawHitIndex>;
  using Mu2eKinKal::Mu2eConfig;
  using CCHandle = art::ValidHandle<CaloClusterCollection>;
  template <class KTRAJ> class KKFit {
    public:
      // fit configuration
      using KKTRK = KKTrack<KTRAJ>;
      using KKTRKPTR = std::unique_ptr<KKTRK>;
      using PKTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      using TCA = KinKal::ClosestApproach<KTRAJ,Line>;
      using KKHIT = KinKal::Measurement<KTRAJ>;
      using KKMAT = KinKal::Material<KTRAJ>;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KKSTRAWXING = KKStrawXing<KTRAJ>;
      using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
      using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
      using KKCALOHIT = KKCaloHit<KTRAJ>;
      using KKCALOHITPTR = std::shared_ptr<KKCALOHIT>;
      using KKCALOHITCOL = std::vector<KKCALOHITPTR>;
      //      using KKPANELHIT = KKPanelHit<KTRAJ>;
      using MEAS = KinKal::Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      // construct from fit configuration objects
      explicit KKFit(Mu2eConfig const& fitconfig);
      // helper functions used to create components of the fit
      void makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, KKBField const& kkbf, StrawMaterial const& smat,
          PKTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
          KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const;
      void addStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, KKBField const& kkbf, StrawMaterial const& smat,
          KKTRK& kktrk, ComboHitCollection const& chcol, KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const;
      void addStraws(Tracker const& tracker, StrawMaterial const& smat, KKTRK& kktrk, KKSTRAWXINGCOL& exings) const;
      bool makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PKTRAJ const& pktraj, KKCALOHITCOL& hits) const;
      void addCaloHit(Calorimeter const& calo, KKTRK& kktrk, CCHandle cchandle, KKCALOHITCOL& hits) const;
      KalSeed createSeed(KKTRK const& kktrk, TrkFitFlag const& seedflag, std::set<double> const& tsave) const;
      TimeRange range(KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) const; // time range from a set of hits and element Xings
      double zTime(PKTRAJ const& trak, double zpos) const; // find the time the trajectory crosses the plane perp to z at the given z position
      bool useCalo() const { return usecalo_; }
      PDGCode::type fitParticle() const { return tpart_;}
      TrkFitDirection fitDirection() const { return tdir_;}
      bool addMaterial() const { return addmat_; }
    private:
      void fillTrackerInfo(Tracker const& tracker) const;
      PDGCode::type tpart_;
      TrkFitDirection tdir_;
      float nullvscale_;
      bool addmat_, usecalo_;
      double caloDt_; // calo time offset; should come from proditions FIXME!
      double caloPosRes_; // calo cluster transverse position resolution; should come from proditions or CaloCluster FIXME!
      double caloPropSpeed_; // effective light propagation speed in a crystal (including reflections).  Should come from prodtions FIXME
      double minCaloEnergy_; // minimum CaloCluster energy
      double maxCaloDt_; // maximum track-calo time difference
      double maxCaloDoca_; // maximum track-calo DOCA
      double tprec_; // TPOCA precision
      StrawHitFlag addsel_, addrej_;
      float maxStrawHitDoca_, maxStrawHitDt_, maxStrawHitChi_, maxStrawDoca_;
      int sbuff_;
      // cached info computed from the tracker; these should be set on construction, but the tracker doesn't exist
      mutable double ymin_, ymax_, umax_; // panel-level info
      mutable double rmin_, rmax_; // plane-level info
      mutable double spitch_;
      mutable bool needstrackerinfo_;
  };

  template <class KTRAJ> KKFit<KTRAJ>::KKFit(Mu2eConfig const& fitconfig) :
    tpart_(static_cast<PDGCode::type>(fitconfig.fitParticle())),
    tdir_(static_cast<TrkFitDirection::FitDirection>(fitconfig.fitDirection())),
    nullvscale_(fitconfig.nullHitVarianceScale()),
    addmat_(fitconfig.addMaterial()),
    usecalo_(fitconfig.useCaloCluster()),
    caloDt_(fitconfig.caloDt()),
    caloPosRes_(fitconfig.caloPosRes()),
    caloPropSpeed_(fitconfig.caloPropSpeed()),
    minCaloEnergy_(fitconfig.minCaloEnergy()),
    maxCaloDt_(fitconfig.maxCaloDt()),
    maxCaloDoca_(fitconfig.maxCaloDoca()),
    tprec_(fitconfig.tpocaPrec()),
    addsel_(fitconfig.addHitSelect()),
    addrej_(fitconfig.addHitReject()),
    maxStrawHitDoca_(fitconfig.maxStrawHitDOCA()),
    maxStrawHitDt_(fitconfig.maxStrawHitDt()),
    maxStrawHitChi_(fitconfig.maxStrawHitChi()),
    maxStrawDoca_(fitconfig.maxStrawDOCA()),
    sbuff_(fitconfig.strawBuffer()),
    needstrackerinfo_(true)
  {
  }

  template <class KTRAJ> void KKFit<KTRAJ>::makeStrawHits(Tracker const& tracker,StrawResponse const& strawresponse,KKBField const& kkbf, StrawMaterial const& smat,
      PKTRAJ const& ptraj, ComboHitCollection const& chcol, StrawHitIndexCollection const& strawHitIdxs,
      KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const {
    // initialize hits as null (no drift).  Drift is turned on when updating
    auto const& sprop = tracker.strawProperties();
    double rstraw = sprop.strawInnerRadius();
    WireHitState whstate(WireHitState::null); // initial wire hit state

    // loop over the individual straw combo hits
    for(auto strawidx : strawHitIdxs) {
      const ComboHit& combohit(chcol.at(strawidx));
      if(combohit.mask().level() != StrawIdMask::uniquestraw){
        throw cet::exception("RECO")<<"mu2e::KKFit: ComboHit error"<< endl;
      }
      const Straw& straw = tracker.getStraw(combohit.strawId());
      auto wline = Mu2eKinKal::hitLine(combohit,straw,strawresponse); // points from the signal to the straw center
      double htime = combohit.correctedTime();
      // use the hit z position to estimate the particle time.
      auto hpos = wline.position3(htime);
      auto ppos = ptraj.position3(htime);
      auto vel = ptraj.velocity(htime);
      double ptime = htime + (hpos.Z()-ppos.Z())/vel.Z();
      CAHint hint(ptime,htime);
      // compute PTCA between the seed trajectory and this straw
      PTCA ptca(ptraj, wline, hint, tprec_ );
      // create the material crossing
      if(addmat_) exings.push_back(std::make_shared<KKSTRAWXING>(ptca,smat,straw.id()));
      // create the hit
      hits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, ptca, whstate, rstraw, combohit, straw, strawidx, strawresponse));
    }
  }

  template <class KTRAJ> bool KKFit<KTRAJ>::makeCaloHit(CCPtr const& cluster, Calorimeter const& calo, PKTRAJ const& ptraj, KKCALOHITCOL& hits) const {
    bool retval(false);
    // move cluster COG into the tracker frame.  COG is at the front face of the disk
    CLHEP::Hep3Vector cog = calo.geomUtil().mu2eToTracker(calo.geomUtil().diskFFToMu2e( cluster->diskID(), cluster->cog3Vector()));
    // project this along the crystal axis to the SIPM, which is at the back.  This is the point the time measurement corresponds to
    VEC3 ffcog(cog);
    double lcrystal = calo.caloInfo().getDouble("crystalZLength");
    VEC3 crystalF2B = VEC3(0.0,0.0,lcrystal); // this should come directly from the calogeometry, TODO
    VEC3 sipmcog = ffcog + crystalF2B;
    // create the Line trajectory from this information: signal goes towards the sipm
    Line caxis(sipmcog,ffcog,cluster->time()+caloDt_,caloPropSpeed_);
    // find the time the seed traj passes the middle of the crystal
    double zt = zTime(ptraj,0.5*(sipmcog.Z()+ffcog.Z()));
    CAHint hint( zt, caxis.t0());
    // compute a preliminary PTCA between the seed trajectory and this straw.
    PTCA ptca(ptraj, caxis, hint, tprec_ );
    // check that this is within tolerance
    if(ptca.usable() && fabs(ptca.doca()) < maxCaloDoca_){
      // check that the sensor position is within the active position of the crystal
      auto poca = ptca.sensorPoca().Vect();
      if((poca.Z()-ffcog.Z()) > -maxCaloDoca_ &&  (poca.Z()-sipmcog.Z()) < maxCaloDoca_) {
        // create the hit
        double tvar = cluster->timeErr()*cluster->timeErr();
        double wvar = caloPosRes_*caloPosRes_;
        hits.push_back(std::make_shared<KKCALOHIT>(cluster,ptca,tvar,wvar));
        retval = true;
      }
    }
    return retval;
  }

  template <class KTRAJ> void KKFit<KTRAJ>::addStrawHits(Tracker const& tracker,StrawResponse const& strawresponse, KKBField const& kkbf, StrawMaterial const& smat,
      KKTRK& kktrk, ComboHitCollection const& chcol,
      KKSTRAWHITCOL& hits, KKSTRAWXINGCOL& exings) const {
    // initialize hits as null (no drift).  Drift is turned on when updating
    auto const& sprop = tracker.strawProperties();
    double rstraw = sprop.strawInnerRadius();
    WireHitState whstate(WireHitState::null); // initial wire hit state;
    auto const& ftraj = kktrk.fitTraj();
    // build the set of existing hits
    std::set<StrawHitIndex> oldhits;
    for(auto const& strawhit : kktrk.strawHits())oldhits.insert(strawhit->strawHitIndex());
    std::set<StrawId> oldstraws;
    for(auto const& strawxing : kktrk.strawXings())oldstraws.insert(strawxing->strawId());
    for(auto const& strawxing : exings)oldstraws.insert(strawxing->strawId());
    for( size_t ich=0; ich < chcol.size();++ich){
      if(oldhits.find(ich)==oldhits.end()){      // make sure this hit wasn't already found
        ComboHit const& strawhit = chcol[ich];
        if(strawhit.flag().hasAllProperties(addsel_) && (!strawhit.flag().hasAnyProperty(addrej_))){
          double zt = zTime(kktrk.fitTraj(),strawhit.pos().Z());
          if(fabs(strawhit.correctedTime()-zt) < maxStrawHitDt_) {      // compare the measured time with the estimate from the fit
            const Straw& straw = tracker.getStraw(strawhit.strawId());
            auto wline = Mu2eKinKal::hitLine(strawhit,straw,strawresponse);
            double psign = wline.direction().Dot(straw.wireDirection());  // wire distance is WRT straw center, in the nominal wire direction
            double htime = wline.t0() - (straw.halfLength()-psign*strawhit.wireDist())/wline.speed();
            CAHint hint(zt,htime);
            // compute PTCA between the trajectory and this straw
            PTCA ptca(ftraj, wline, hint, tprec_ );
            if(fabs(ptca.doca()) < maxStrawHitDoca_){ // add test of chi TODO
              hits.push_back(std::make_shared<KKSTRAWHIT>(kkbf, ptca, whstate, rstraw, strawhit, straw, ich, strawresponse));
              if(addmat_ && oldstraws.find(straw.id()) == oldstraws.end())exings.push_back(std::make_shared<KKSTRAWXING>(ptca,smat,straw.id()));
            }
          }
        }
      }
    }
  }

  template <class KTRAJ> void KKFit<KTRAJ>::addStraws(Tracker const& tracker, StrawMaterial const& smat, KKTRK& kktrk, KKSTRAWXINGCOL& exings) const {
    // this algorithm assumes the track never hits the same straw twice.  That could be violated by reflecting tracks, and could be addressed
    // by including the time of the Xing as part of its identity.  That would slow things down so it remains to be proven it's a problem  TODO
    // build the set of existing straws
    if(addmat_){
      auto const& ftraj = kktrk.fitTraj();
      // pre-compute some tracker info if needed
      if(needstrackerinfo_)fillTrackerInfo(tracker);
      std::set<StrawId> oldstraws;
      for(auto const& strawxing : kktrk.strawXings())oldstraws.insert(strawxing->strawId());
      for(auto const& strawxing : exings)oldstraws.insert(strawxing->strawId());
      // test for the track going through a panel.  Start with planes
      for(auto const& plane : tracker.planes()){
        if(tracker.planeExists(plane.id())) {
          double plz = plane.origin().z();
          // find the track position in at this plane's central z.
          double zt = zTime(ftraj,plz);
          auto plpos = ftraj.position3(zt);
          // rough check on the point radius
          double rho = plpos.Rho();
          if(rho > rmin_ && rho < rmax_){
            // loop over panels in this plane
            for(auto panel_p : plane.panels()){
              auto const& panel = *panel_p;
              // linearly correct position for the track direction due to difference in panel-plane Z position
              auto tdir = ftraj.direction(zt);
              double dz = panel.origin().z()-plz;
              auto papos = plpos + (dz/tdir.Z())*tdir;
              // convert this position into panel coordinates
              CLHEP::Hep3Vector cpos(papos.X(),papos.Y(),papos.Z()); // clumsy translation
              auto pposv = panel.dsToPanel()*cpos;
              // translate the y position into a rough straw number
              int istraw = static_cast<int>(rint( (pposv.y()-ymin_)*spitch_));
              // require this be within the (integral) straw buffer
              if(istraw >= -sbuff_ && istraw < static_cast<int>(panel.nStraws()) + sbuff_ ){
                unsigned istrmin = static_cast<unsigned>(std::max(istraw-sbuff_,0));
                // largest straw is the innermost; use that to test length
                if(fabs(pposv.x()) < panel.getStraw(istrmin).halfLength() ) {
                  unsigned istrmax = static_cast<unsigned>(std::min(istraw+sbuff_,static_cast<int>(panel.nStraws())-1));
                  // loop over straws
                  for(unsigned istr = istrmin; istr <= istrmax; ++istr){
                    auto const& straw = panel.getStraw(istr);
                    // add strawExists test TODO
                    // make sure we haven't already seen this straw
                    if(oldstraws.find(straw.id()) == oldstraws.end()){
                      KinKal::VEC3 vp0(straw.wireEnd(StrawEnd::cal));
                      KinKal::VEC3 vp1(straw.wireEnd(StrawEnd::hv));
                      KinKal::VEC3 smid = 0.5*(vp0+vp1);
                      KinKal::Line wline(vp0,vp1,zt,CLHEP::c_light); // time is irrelevant: use speed of light as sprop
                      CAHint hint(zt,zt);
                      // compute PTCA between the trajectory and this straw
                      PTCA ptca(ftraj, wline, hint, tprec_ );
                      double du = (ptca.sensorPoca().Vect()-smid).R();
                      if(fabs(ptca.doca()) < maxStrawDoca_ && du < straw.halfLength()){ // add test of chi TODO
                        exings.push_back(std::make_shared<KKSTRAWXING>(ptca,smat,straw.id()));
                      }
                    } // not existing straw cut
                  } // straws loop
                } // straw length cut
              } // straw index cut
            } // panels loop
          } // radius cut
        } // plane exists cut
      }  // planes loop
    } // adding material
  } // end function

  template <class KTRAJ> void KKFit<KTRAJ>::addCaloHit(Calorimeter const& calo, KKTRK& kktrk, CCHandle cchandle, KKCALOHITCOL& hits) const {
    double crystalLength = calo.caloInfo().getDouble("crystalZLength");
    auto const& ftraj = kktrk.fitTraj();
    bool added(false);
    auto cccol = cchandle.product();
    unsigned idisk=0;
    // loop over disks, or until we find a matching cluster
    while(!added && idisk < 2){
      auto ffpos = calo.geomUtil().mu2eToTracker(calo.disk(idisk).geomInfo().frontFaceCenter());
      double rmin = calo.disk(idisk).geomInfo().innerEnvelopeR();
      double rmax = calo.disk(idisk).geomInfo().outerEnvelopeR();
      // test at both faces; if the track is in the right area, test the clusters on this disk
      bool test(false);
      for(int iface=0; iface<2; ++iface){
        double zt = zTime(ftraj,ffpos.z()+iface*crystalLength);
        auto tpos = ftraj.position3(zt);
        double rho = tpos.Rho();
        if ( rho > rmin && rho < rmax ){
          test=true;
          break;
        }
      }
      if(test){
        for(size_t icc=0;icc < cccol->size(); ++icc){
          auto const& cc = (*cccol)[icc];
          if (static_cast<unsigned>(cc.diskID()) == idisk && cc.energyDep() > minCaloEnergy_){
            auto ccPtr = art::Ptr<CaloCluster>(cchandle,icc);
            added = makeCaloHit(ccPtr,calo,ftraj,hits);
            if(added)break;
          }
        }
      }
      ++idisk;
    }
  }

  template <class KTRAJ> void KKFit<KTRAJ>::fillTrackerInfo(Tracker const& tracker) const {
    // pre-compute some info derived from the tracker
    double strawradius = tracker.strawOuterRadius();
    auto const& frontplane = tracker.planes().front();
    auto const& firstpanel = frontplane.getPanel(0);
    auto const& innerstraw = firstpanel.getStraw(0);
    auto const& outerstraw = firstpanel.getStraw(StrawId::_nstraws-1);
    // compute limits: add some buffer for the finite size of the straw
    auto DStoP = firstpanel.dsToPanel();
    auto innerstraw_origin = DStoP*innerstraw.origin();
    auto outerstraw_origin = DStoP*outerstraw.origin();
    ymin_ = innerstraw_origin.y();
    ymax_ = outerstraw_origin.y();
    umax_ = innerstraw.halfLength() + strawradius; // longest possible straw
    // plane-level variables: these add some buffer
    rmin_ = innerstraw_origin.y() - sbuff_*strawradius;
    rmax_ = outerstraw.wireEnd(StrawEnd::cal).mag() + sbuff_*strawradius;
    spitch_ = (StrawId::_nstraws-1)/(ymax_-ymin_);
    needstrackerinfo_= false;
  }

  template <class KTRAJ> double KKFit<KTRAJ>::zTime(PKTRAJ const& ptraj, double zpos) const {
    if(ptraj.range().null())throw cet::exception("RECO")<<"mu2e::KKFit: PTraj range error in zTime"<< endl;
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

  template <class KTRAJ> TimeRange KKFit<KTRAJ>::range(KKSTRAWHITCOL const& strawhits, KKCALOHITCOL const& calohits, KKSTRAWXINGCOL const& strawxings) const{
    double tmin = std::numeric_limits<float>::max();
    double tmax = std::numeric_limits<float>::min();
    for( auto const& strawhit : strawhits) {
      tmin = std::min(tmin,strawhit->time());
      tmax = std::max(tmax,strawhit->time());
    }
    for( auto const& calohit : calohits) {
      tmin = std::min(tmin,calohit->time());
      tmax = std::max(tmax,calohit->time());
    }
    for( auto const& strawxing : strawxings) {
      tmin = std::min(tmin,strawxing->crossingTime());
      tmax = std::max(tmax,strawxing->crossingTime());
    }
    return TimeRange(tmin,tmax);
  }

  template <class KTRAJ> KalSeed KKFit<KTRAJ>::createSeed(KKTRK const& kktrk, TrkFitFlag const& seedflag, std::set<double> const& savetimes) const {
    TrkFitFlag fflag(seedflag);  // initialize the flag with the seed fit flag
    if(kktrk.fitStatus().usable())fflag.merge(TrkFitFlag::kalmanOK);
    if(kktrk.fitStatus().status_ == Status::converged) fflag.merge(TrkFitFlag::kalmanConverged);
    if(addmat_)fflag.merge(TrkFitFlag::MatCorr);
    if(kktrk.config().bfcorr_ )fflag.merge(TrkFitFlag::BFCorr2);
    // explicit T0 is needed for backwards-compatibility; sample from the appropriate trajectory piece
    auto const& fittraj = kktrk.fitTraj();
    double tz0 = zTime(fittraj,0.0);
    auto const& t0piece = fittraj.nearestPiece(tz0);
    HitT0 t0(t0piece.paramVal(KTRAJ::t0_), sqrt(t0piece.paramVar(KTRAJ::t0_)));
    // create the shell for the output.  Note the (obsolete) flight length is given as t0
    KalSeed fseed(tpart_,tdir_,fflag,t0.t0());
    auto const& fstatus = kktrk.fitStatus();
    fseed._chisq = fstatus.chisq_.chisq();
    fseed._fitcon = fstatus.chisq_.probability();
    fseed._nseg = fittraj.pieces().size();
    // loop over track components and store them
    fseed._hits.reserve(kktrk.strawHits().size());
    for(auto const& strawhit : kktrk.strawHits() ) {
      auto const& chit = strawhit->hit();
      StrawHitFlag hflag = chit.flag();
      if(strawhit->active())hflag.merge(StrawHitFlag::active);
      auto const& ca = strawhit->closestApproach();
      TrkStrawHitSeed seedhit(strawhit->strawHitIndex(),
          HitT0(ca.particleToca(), sqrt(ca.tocaVar())),
          static_cast<float>(ca.particleToca()), static_cast<float>(ca.sensorToca()),
          static_cast<float>(strawhit->timeResidual().value()), // drift radius isn't a clear concept in KinKal, use time residual instead
          static_cast<float>(strawhit->time()),
          static_cast<float>(ca.doca()),
          strawhit->hitState().state_,
          static_cast<float>(sqrt(strawhit->timeResidual().variance())), // substitute for drift radius error
          hflag, chit);
      fseed._hits.push_back(seedhit);
    }
    if(kktrk.caloHits().size() > 0){
      auto const& calohit = kktrk.caloHits().front(); // for now take the front: not sure if there will ever be >1 TODO
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
    }
    fseed._straws.reserve(kktrk.strawXings().size());
    for(auto const& sxing : kktrk.strawXings()) {
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      // compute energy loss
      sxing->materialEffects(kktrk.refTraj(),KinKal::TimeDir::forwards, dmom, momvar);
      VEC3 dm(dmom[0],dmom[1],dmom[2]);
      // find material crossing properties
      double gaspath(0.0);
      double radfrac(0.0);
      for(auto const& matxing : sxing->matXings()){
        if(matxing.dmat_.name().compare(5,3,"gas",3)==0)gaspath = matxing.plen_;
        radfrac += matxing.dmat_.radiationFraction(matxing.plen_);
        auto const& tpocad = sxing->closestApproach();
        fseed._straws.emplace_back(sxing->strawId(),
            tpocad.doca(),
            tpocad.particleToca(),
            tpocad.sensorToca(),
            gaspath,
            radfrac,
            dm.R(),
            sxing->active() );
      }
    }
    // sample the fit at the requested times and save those segments.  Uniqueness needs to be checked in the calling function
    fseed._segments.reserve(savetimes.size());
    for(auto time : savetimes) fseed._segments.emplace_back(fittraj.nearestPiece(time),time);
    return fseed;
  }

}
#endif
