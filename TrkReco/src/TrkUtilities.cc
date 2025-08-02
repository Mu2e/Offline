//
//  Collection of tools useful for dealing with various helices
//  Original Author Dave Brown (LBNL) 26 Aug. 2016
//
// Mu2e
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkStraw.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
// BTrk
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "Offline/BTrkData/inc/TrkCaloHit.hh"
#include "Offline/Mu2eBTrk/inc/DetStrawElem.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
// KinKal
#include "KinKal/Trajectory/CentralHelix.hh"
// CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
//C++
#include <cmath>
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e {
  namespace TrkUtilities {

    bool RobustHelix2Traj (RobustHelix const& helix, HepVector& hpvec, float amsign) {
      bool retval(false);
      // compare the input with this configuration's helicity: these must be the same
      // radius = 0 or lambda=0 are degenerate cases that this representation can't handle
      if(helix.radius() > 0.0 && helix.lambda() != 0.0 && hpvec.num_row() == HelixTraj::NHLXPRM) {
        // radius and omega have inverse magnitude, omega is signed by the angular momentum
        hpvec[HelixTraj::omegaIndex] = amsign/helix.radius();
        // phi0 is the azimuthal angle of the particle velocity vector at the point
        // of closest approach to the origin.  It's sign also depends on the angular
        // momentum.  To translate from the center, we need to reverse coordinates
        hpvec[HelixTraj::phi0Index] = atan2(-amsign*helix.centerx(),amsign*helix.centery());
        // d0 describes the distance to the origin at closest approach.
        // It is signed by the particle angular momentum WRT the origin.
        // The Helix fit radial bias is anti-correlated with d0; correct for it here.
        hpvec[HelixTraj::d0Index] = amsign*(helix.rcent() - helix.radius());
        // the dip angle is measured WRT the perpendicular, signed by the z component of linear momentum
        hpvec[HelixTraj::tanDipIndex] = amsign*helix.lambda()/helix.radius();
        // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
        // resolve the loop ambiguity such that the POCA is closest to z=0.
        double refphi = helix.fz0()+amsign*M_PI_2;
        double phi = hpvec[HelixTraj::phi0Index];
        double dphi = Angles::deltaPhi(phi,refphi);
        // choose z0 (which loop) so that f=0 is as close to z=0 as possible
        hpvec[HelixTraj::z0Index] = dphi*hpvec[HelixTraj::tanDipIndex]/hpvec[HelixTraj::omegaIndex];
        retval = true;
      }
      return retval;
    }

    void RobustHelixFromMom(Hep3Vector const& pos, Hep3Vector const& mom, double charge, double Bz, RobustHelix& helix){
      double momToRad = 1000.0/(charge*Bz*CLHEP::c_light);
      // compute some simple useful parameters
      double pt = mom.perp();
      // transverse radius of the helix
      helix._radius = fabs(pt*momToRad);
      //longitudinal wavelength; sign convention goes with angular rotation
      helix._lambda = -mom.z()*momToRad;
      // circle center
      Hep3Vector center = Hep3Vector(pos.x() + mom.y()*momToRad, pos.y() - mom.x()*momToRad, 0.0);
      helix._rcent = center.perp();
      helix._fcent = center.phi();
      // phi at z=0
      double phi = (pos - center).phi() - pos.z()/helix.lambda();
      // reset to be close to 0
      Angles::deltaPhi(phi);
      helix._fz0 = phi;
    }

    // legacy function
    void fillSegment(HelixTraj const& htraj, double locflt, double globflt, TrkT0 t0, double mass, int charge, BField const& bfield, KalSegment& kseg) {
      // compute the kinematics; this is external to htraj
      double radToMom = charge*bfield.bFieldNominal()*CLHEP::c_light/1000.0;
      double mom = fabs(radToMom/(htraj.omega()*htraj.cosDip()));
      double energy = sqrt(mom*mom + mass*mass);
      double v = CLHEP::c_light*mom/energy;
      double vz = v*htraj.sinDip();
      // translate BTrk t0 to CentralHelix t0 (different convention)
      double ct0 = t0.t0() + htraj.z0()/vz;
      // translate htraj (3D) flight range to time ranges
      double tmin = ct0 + htraj.lowRange()/v;
      double tmax = ct0 + htraj.hiRange()/v;
      double tref = ct0 + locflt/v;
      // conver the helix content to a CentralHelix.  Note the t0 value supplied is in the BTrk convention (time at z=0).
      KinKal::DVEC chpars;
      KinKal::DMAT cov;
      for(unsigned ipar=0; ipar<5; ipar++){
        chpars(ipar) = htraj.parameters()->parameter()[ipar];
        for(unsigned jpar=0; jpar<5; jpar++){
          cov(ipar,jpar) = htraj.parameters()->covariance().fast(ipar+1,jpar+1);
        }
      }
      // insert t0 by hand
      chpars(KinKal::CentralHelix::t0_) = ct0;
      cov(KinKal::CentralHelix::t0_,KinKal::CentralHelix::t0_) = t0.t0Err()*t0.t0Err();
      KinKal::Parameters params(chpars,cov);
      // create the CentralHelix from these
      KinKal::CentralHelix chelix(params,mass,charge,bfield.bFieldNominal(),KinKal::TimeRange(tmin,tmax));
      kseg = KalSegment(chelix, tref, globflt-locflt);
    }

    void fillStraws(const KalRep* krep, std::vector<TrkStraw>& tstraws) {
      tstraws.clear();
// disabled: not compatible with KinKal

      // get material sites from the KalRep
//      for(auto isite : krep->siteList()){
//        if(isite->kalMaterial() != 0) {
//          const KalMaterial* kmat = isite->kalMaterial();
//          const DetStrawElem* detstraw = dynamic_cast<const DetStrawElem*>(kmat->detElem());
//          if(detstraw != 0){
//            // found a straw: create a TrkStraw object from it
//            // i must recompute POCA since the KalMaterial doesn't cache the hit flight FIXME!
//            TrkPoca poca(krep->traj(),kmat->detIntersection().pathlen,*detstraw->wireTraj(),0);
//            TrkStraw tstraw(detstraw->straw()->id(),
//                kmat->detIntersection().dist, //poca.doca(),
//                kmat->detIntersection().pathlen, // poca.flt1(),
//                poca.flt2(),  // not stored in KalMaterial, FIXME!
//                kmat->detIntersection().pathLength(),
//                detstraw->radiationFraction(kmat->detIntersection()),
//                kmat->momFraction(),
//                isite->isActive() );
//            tstraws.push_back(tstraw);
//          }
//        }
//      }
    }

    void fillStrawHitSeeds(const KalRep* krep,ComboHitCollection const& chits, std::vector<TrkStrawHitSeed>& hitseeds) {
      // extract the TkrStrawHits from the KalRep
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      // loop over the TrkStrawHits and convert them
      for(auto tsh : tshv ) {
        // find the associated ComboHit
        auto const& chit = chits.at(tsh->index());
        // set the flag according to the status of this hit
        StrawHitFlag hflag = chit.flag();
        if(tsh->isActive())hflag.merge(StrawHitFlag::active);
        if(tsh->poca().status().success())hflag.merge(StrawHitFlag::doca);
        int state = tsh->ambig();
        if(!tsh->isActive())state = -2;
        CLHEP::Hep3Vector hpos = tsh->hitTraj()->position(tsh->hitLen());
        float upos = (hpos-tsh->straw().wirePosition(0.0)).dot(tsh->straw().wireDirection());
        float dt = tsh->signalTime() - tsh->driftTime()-tsh->hitT0()._t0;
        hitseeds.emplace_back(tsh->index(),
            tsh->hitT0(), tsh->fltLen(), tsh->hitLen(),
            tsh->driftRadius(), tsh->signalTime(), upos, dt,
            tsh->poca().doca(), state, tsh->driftRadiusErr(), hflag, chit);
      }
    }


    void fillCaloHitSeed(const TrkCaloHit* tch, CLHEP::Hep3Vector const& tmom, TrkCaloHitSeed& caloseed) {
      // set the flag according to the status of this hit
      StrawHitFlag hflag;
      if(tch->isActive())hflag.merge(StrawHitFlag::active);
      if(tch->poca().status().success())hflag.merge(StrawHitFlag::doca);
      Hep3Vector hpos;
      tch->hitPosition(hpos);
      caloseed = TrkCaloHitSeed(tch->hitT0(), tch->fltLen(), tch->hitLen(),
          tch->poca().doca(), tch->hitErr(), tch->time() + tch->timeOffset(), tch->timeErr(),
          XYZVectorF(hpos),
          XYZVectorF(tmom),
          hflag);
    }
    // DNB: the timeOffset() should NOT be added to time(), it is a double correction.
    // I'm leaving for now as the production was run with this error FIXME!

    // compute the overlap between 2 clusters
    double overlap(SHIV const& shiv1, SHIV const& shiv2) {
      double over(0.0);
      double norm = std::min(shiv1.size(),shiv2.size());
      // count the overlapping hits
      for(auto h1 : shiv1){
        for(auto h2 : shiv2){
          if(h1 == h2){
            over +=1.0;
            break;
          }
        }
      }
      return over/norm;
    }

    double overlap(TimeCluster const& tc1, TimeCluster const& tc2) {
      double hover = overlap(tc1._strawHitIdxs,tc2._strawHitIdxs);
      double norm = std::min(tc1.hits().size(),tc2.hits().size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(tc1._caloCluster.isNonnull() && tc2._caloCluster.isNonnull()) {
        if(tc1._caloCluster == tc2._caloCluster)
          over += norm;
        norm *= 2;
      }
      return over/norm;
    }

    double overlap(KalSeed const& ks1, KalSeed const& ks2) {
      // translate hit info into a simple index array.  Only count active hits
      SHIV shiv1, shiv2;
      for(auto tshs : ks1.hits()){
        if(tshs.flag().hasAllProperties(StrawHitFlag::active))
          shiv1.push_back(tshs.index());
      }
      for(auto tshs : ks2.hits()){
        if(tshs.flag().hasAllProperties(StrawHitFlag::active))
          shiv2.push_back(tshs.index());
      }
      double hover = overlap(shiv1,shiv2);
      double norm = std::min(shiv1.size(),shiv2.size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(ks1.caloCluster().isNonnull() && ks2.caloCluster().isNonnull()) {
        if(ks1.caloCluster() == ks2.caloCluster())
          over += norm;
        norm *= 2;
      }
      return over/norm;
    }

    double overlap(KalSeed const& ks, HelixSeed const& hs){
      SHIV shiv1, shiv2;
      for(auto tshs : ks.hits()){
        if(tshs.flag().hasAllProperties(StrawHitFlag::active))
          shiv1.push_back(tshs.index());
      }
      for(auto hh : hs.hits()){
        // exclude outliers
        if(!hh.flag().hasAnyProperty(StrawHitFlag::outlier))
          shiv2.push_back(hh.index());
      }
      double hover = overlap(shiv1,shiv2);
      double norm = std::min(shiv1.size(),shiv2.size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(ks.caloCluster().isNonnull() && hs.caloCluster().isNonnull()) {
        if(ks.caloCluster() == hs.caloCluster())
          over += norm;
        norm *= 2;
      }
      return over/norm;
    }

    double overlap(HelixSeed const& hs,TimeCluster const& tc) {
      SHIV shiv;
      for(auto hh : hs.hits()){
        // exclude outliers
        if(!hh.flag().hasAnyProperty(StrawHitFlag::outlier))
          shiv.push_back(hh.index());
      }
      double hover = overlap(shiv,tc.hits());
      double norm = std::min(shiv.size(),tc.hits().size());
      double over = norm*hover;
      // add in CaloCluster; count is as much as all the hits
      if(tc.caloCluster().isNonnull() && hs.caloCluster().isNonnull()) {
        if(tc.caloCluster() == hs.caloCluster())
          over += norm;
        norm *= 2;
      }
      return over/norm;
    }

    // this function belongs in TrkDifTraj, FIXME!!!!
    // double zFlight(TrkDifPieceTraj const& ptraj, double pz) {
    //   // get the helix at the middle of the track
    //   double loclen;
    //   double fltlen(0.0);
    //   const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(ptraj.localTrajectory(fltlen,loclen));
    //   // Iterate
    //   const HelixTraj* oldtraj;
    //   unsigned iter(0);
    //   do {
    //  // remember old traj
    //  oldtraj = htraj;
    //  // correct the global fltlen for this difference in local trajectory fltlen at this Z position
    //  fltlen += (htraj->zFlight(pz)-loclen);
    //  htraj = dynamic_cast<const HelixTraj*>(ptraj.localTrajectory(fltlen,loclen));
    //   } while(oldtraj != htraj && iter++<10);
    //   return fltlen;
    // }

    void countHits(const std::vector<TrkStrawHitSeed>& hits, unsigned& nhits, unsigned& nactive, unsigned& ndouble, unsigned& ndactive, unsigned& nnullambig) {
      nhits = 0; nactive = 0; ndouble = 0; ndactive = 0; nnullambig = 0;
      static StrawHitFlag active(StrawHitFlag::active);
      for (std::vector<TrkStrawHitSeed>::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit) {
        ++nhits;
        if (ihit->flag().hasAllProperties(active)) {
          ++nactive;
          if (ihit->ambig()==0) {
            ++nnullambig;
          }
        }
        /*    if (ihit->nStrawHits()>=2) {
              ++ndactive;
              }
              */
        //    std::cout << "AE: ihit->nStrawHits() = " << ihit->nStrawHits() << std::endl;
        const auto& jhit = ihit+1;
        const auto& hhit = ihit-1;
        if( (jhit != hits.end() &&
              jhit->flag().hasAllProperties(active) &&
              jhit->strawId().getPlane() == ihit->strawId().getPlane() &&
              jhit->strawId().getPanel() == ihit->strawId().getPanel() ) ||
            (hhit >= hits.begin() &&
             hhit->flag().hasAllProperties(active) &&
             hhit->strawId().getPlane() == ihit->strawId().getPlane() &&
             hhit->strawId().getPanel() == ihit->strawId().getPanel() )
          ) {
          ++ndouble;
          if (ihit->flag().hasAllProperties(StrawHitFlag::active)) {
            ++ndactive;
          }
        }
      }
      //      std::cout << "AE: ndactive hits = " << ndactive << std::endl;
    }

    double chisqConsistency(const KalRep* krep) {
      return ChisqConsistency(krep->chisq(),krep->nDof()-1).significanceLevel();
    }

    unsigned countBends(const KalRep* krep) {
      unsigned nbend(0);
      for(auto isite : krep->siteList()){
        if(isite->kalBend() != 0) ++nbend;
      }
      return nbend;
    }

    const TrkCaloHit* findTrkCaloHit(const KalRep* krep){
      const TrkCaloHit* tch(0);
      for(auto ith=krep->hitVector().begin(); ith!=krep->hitVector().end(); ++ith){
        const TrkCaloHit* tsh = dynamic_cast<const TrkCaloHit*>(*ith);
        if(tsh != 0) {
          tch = tsh;
          break;
        }
      }
      return tch;
    }
    double energy(double mass, double momentum) {  return sqrt(momentum*momentum + mass*mass); }
    double beta(double mass, double momentum) { return fabs(momentum)/energy(mass,momentum); }
    double betagamma(double mass, double momentum) { return fabs(momentum)/mass; }
    double gamma(double mass, double momentum) { return energy(mass,momentum)/mass; }
  } // TrkUtilities
}// mu2e
