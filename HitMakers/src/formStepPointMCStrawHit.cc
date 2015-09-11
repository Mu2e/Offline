//
// Free function to form StepPointMCStrawHit's
//
// $Id: formStepPointMCStrawHit.cc,v 1.2 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author KLG based on Rob's MakeStrawHit_module
//

// Mu2e includes
#include "HitMakers/inc/formStepPointMCStrawHit.hh"
#include "GeneralUtilities/inc/LinePointPCA.hh"
#include "GeneralUtilities/inc/TwoLinePCA.hh"

namespace mu2e {

  std::unique_ptr<StepPointMCStrawHit> formStepPointMCStrawHit(
                                               art::Ptr<StepPointMC> const spmcp,
                                               StrawIndex const & straw_id,
                                               double _minimumLength,
                                               bool   _enableFlightTimeCorrection,
                                               MassCache & cache,
                                               CLHEP::RandGaussQ & _gaussian,
                                               Tracker const & tracker,
                                               ConditionsHandle<TrackerCalibrations> const & trackerCalibrations
                                               ) {

    StrawIndex spmcStraw_id      = spmcp->strawIndex();
    Straw const&  straw          = tracker.getStraw(spmcStraw_id);
    CLHEP::Hep3Vector const& mid = straw.getMidPoint();
    CLHEP::Hep3Vector const& w   = straw.getDirection();
    double strawHalfLength       = straw.getHalfLength();
    double signalVelocity        = trackerCalibrations->SignalVelocity(spmcStraw_id);

    CLHEP::Hep3Vector  const& pos = spmcp->position();
    CLHEP::Hep3Vector  const& mom = spmcp->momentum();
    double length  = spmcp->stepLength();
    // FIXME we calculate the crosstalk using the energy..., should we convert to amplitude here?
    // or postpone the adjustment?
    double edep = (straw_id == spmcStraw_id) ?
      spmcp->totalEDep() :
      spmcp->totalEDep()*trackerCalibrations->CrossTalk(straw_id,spmcStraw_id);
    double hitTime = spmcp->time();

    // Calculate the drift distance and the point on the wire at the
    // end of the dca vector.

    double hit_dca;
    CLHEP::Hep3Vector hit_pca;

    // Length along the step from the start to the dca.
    double hit_s(0.);

    if( length < _minimumLength ) {

      // If step length is very small, consider it a point

      LinePointPCA pca(mid, w, pos);
      hit_dca = pca.dca();
      hit_pca = pca.pca();
      hit_s   = 0.;

    } else {

      // Step is not a point. Calculate the distance between two lines.

      TwoLinePCA pca( mid, w, pos, mom);

      if ( pca.s2() >=0 && pca.s2() <= length ){

        // If the point of closest approach is within the step and wire - thats it.
        hit_dca = pca.dca();
        hit_pca = pca.point1();
        hit_s   = pca.s2();

      } else {


        // The point of closest approach is not within the step. In this case
        // the closes distance should be calculated from the ends

        LinePointPCA pca1(mid, w, pos);
        LinePointPCA pca2(mid, w, pos+mom.unit()*length);
        if( pca1.dca() < pca2.dca() ) {
          hit_dca = pca1.dca();
          hit_pca = pca1.pca();
          hit_s   = 0.;

        } else {
          hit_dca = pca2.dca();
          hit_pca = pca2.pca();
          hit_s   = length;
        }

      }

    } // drift distance calculation

    // Flight time of particle from start of step to the DCA.
    double flightTime = 0.;
    if ( _enableFlightTimeCorrection ){
      double mass = cache.mass( spmcp->simParticle()->pdgId() );
      double p    = mom.mag();
      double e    = sqrt( p*p + mass*mass);
      double beta = p/e;
      flightTime  = ( beta > 0 ) ? hit_s/beta/CLHEP::c_light : 0.;
    }

    // Calculate signal time. It is Geant4 time + signal propagation time
    // t1 is signal time at positive end (along w vector),
    // t2 - at negative end (opposite to w vector)

    D2T d2t;
    trackerCalibrations->DistanceToTime(spmcStraw_id,hit_dca,mom,d2t);
    // smear the time to account for dispersion and measurement error.  Truncate at 0
    double driftTime = std::max(0.0,d2t._tdrift + _gaussian.fire(0.,d2t._tdrifterr));
    double distanceToMiddle = (hit_pca-mid).dot(w);

    // The convention is that the principal time measurement (t1)
    // corresponds to a measurement at the end of the wire as signed
    // by the wire direction vector. t2 is at the near end.

    double hit_t1 = hitTime+flightTime+driftTime+
      (strawHalfLength-distanceToMiddle)/signalVelocity;
    double hit_t2 = hit_t1 + 2.*distanceToMiddle/signalVelocity;

    return std::unique_ptr<mu2e::StepPointMCStrawHit>(
           new StepPointMCStrawHit(spmcp,
                                   edep,
                                   hit_dca,
                                   d2t._tdrift,
                                   driftTime,
                                   distanceToMiddle,
                                   hit_t1,
                                   hit_t2));

  } //end formStepPointMCStrawHit

}
