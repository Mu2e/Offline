// Author : S Middleton
// Date : 2019
// Purpose : Stores functions associted with Cosmic Drift Fit

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "TMatrixD.h"

//Tracker Drift Conditions:
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/StrawPhysics.hh"
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"

using namespace mu2e;

    double DriftFitUtils::GetTestDOCA(ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker) {

        Straw const& straw = tracker->getStraw(chit.strawId());

        XYZVectorF track_position(a0,b0,0);
        XYZVectorF track_direction(a1,b1,1);

        const CLHEP::Hep3Vector& spos = straw.getMidPoint();
        const CLHEP::Hep3Vector& sdir = straw.getDirection();

        XYZVectorF wire_position = XYZVectorF(spos);
        XYZVectorF wire_direction= XYZVectorF(sdir);

        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
                track_direction,
                wire_position,
                wire_direction,
                1.e-8);

        double dca;
        dca = PCA.dca();

        return dca;
}

    double DriftFitUtils::GetTestDOCA(ComboHit const& chit, XYZVectorF const& track_position, XYZVectorF const&  track_direction, const Tracker* tracker) {

        Straw const& straw = tracker->getStraw(chit.strawId());
        const CLHEP::Hep3Vector& spos = straw.getMidPoint();
        const CLHEP::Hep3Vector& sdir = straw.getDirection();

        XYZVectorF wire_position = XYZVectorF(spos);
        XYZVectorF wire_direction= XYZVectorF(sdir);

        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
                track_direction,
                wire_position,
                wire_direction,
                1.e-8);

        double dca;
        dca = PCA.dca();

        return dca;
}

     double DriftFitUtils::GetRPerp(StrawResponse const& _srep, ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker){
        StrawId const& straw_id = chit.strawId();
        Straw const& straw = tracker->getStraw(straw_id);
        XYZVectorF track_pos(a0,b0,0);
        XYZVectorF track_dir(a1,b1,1);

         Hep3Vector td(a1, b1, 1);
        td = td.unit();
        Hep3Vector rperp = td - (td.dot(straw.getDirection()) * straw.getDirection());
        auto const& straw_mp = straw.getMidPoint();
        auto const& wire_dir = straw.getDirection().unit();
        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_pos, track_dir, XYZVectorF(straw_mp),
                                        XYZVectorF(wire_dir), 1.e-8);

        double phi = rperp.theta();
        double drift_distance =
        _srep.driftTimeToDistance(straw_id, chit.driftTime(),
                                  phi);
        double residual = (PCA.LRambig() * PCA.dca()) - drift_distance;
        return residual;
}

    double DriftFitUtils::GetDriftDistance(StrawResponse const& _srep, ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker){
        StrawId const& straw_id = chit.strawId();
        Straw const& straw = tracker->getStraw(straw_id);
        XYZVectorF track_pos(a0,b0,0);
        XYZVectorF track_dir(a1,b1,1);

         Hep3Vector td(a1, b1, 1);
        td = td.unit();
        Hep3Vector rperp = td - (td.dot(straw.getDirection()) * straw.getDirection());

        double phi = rperp.theta();
        double drift_distance = _srep.driftTimeToDistance(straw_id, chit.driftTime(),
                                  phi);
        return drift_distance;
}


    int DriftFitUtils::GetAmbig(ComboHit const& chit, XYZVectorF const& track_position, XYZVectorF const&  track_direction, const Tracker* tracker) {
        Straw const& straw = tracker->getStraw(chit.strawId());

        const CLHEP::Hep3Vector& spos = straw.getMidPoint();
        const CLHEP::Hep3Vector& sdir = straw.getDirection();

        XYZVectorF wire_position = XYZVectorF(spos);
        XYZVectorF wire_direction= XYZVectorF(sdir);

        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
                track_direction,
                wire_position,
                wire_direction,
                1.e-8);

        double ambig = PCA.LRambig();

              return ambig;

}

    int DriftFitUtils::GetAmbig(ComboHit const& chit, double a0, double a1, double b0, double b1,  const Tracker* tracker) {
        Straw const& straw = tracker->getStraw(chit.strawId());

        XYZVectorF track_position(a0,b0,0);
        XYZVectorF track_direction(a1,b1,1);

        const CLHEP::Hep3Vector& spos = straw.getMidPoint();
        const CLHEP::Hep3Vector& sdir = straw.getDirection();

        XYZVectorF wire_position = XYZVectorF(spos);
        XYZVectorF wire_direction= XYZVectorF(sdir);

        TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
                track_direction,
                wire_position,
                wire_direction,
                1.e-8);

        double ambig = PCA.LRambig();

              return ambig;
}

    double DriftFitUtils::GetPropVelocity(StrawResponse const& rep, ComboHit const& chit){
           double vprop = 2.0*rep.halfPropV(chit.strawId(),1000.0*chit.energyDep());
        return vprop;
}

    double DriftFitUtils::GetPropTime(ComboHit const& chit,  double vprop, const Tracker* tracker) {

        Straw const& straw = tracker->getStraw(chit.strawId());

        double tprop = 0.;
        switch (chit.earlyEnd()) {
        case StrawEnd::cal:
          tprop = (straw.halfLength()+chit.wireDist())/vprop;
          break;
        case StrawEnd::hv:
          tprop = (straw.halfLength()-chit.wireDist())/vprop;
          break;
        }

        return tprop;
    }


   double DriftFitUtils::TimeResidualTrans( double doca){
        double drift_time= doca/0.065;
        return drift_time;
   }

    double DriftFitUtils::TimeResidualLong(double doca, StrawResponse const& srep,  double t0, ComboHit const& chit,  const Tracker* tracker){
        double _vprop = 2.0*srep.halfPropV(chit.strawId(),1000.0*chit.energyDep());
        double propagation_time = GetPropTime(chit, _vprop, tracker);
        return propagation_time;
}

    double DriftFitUtils::TimeResidual(double doca, StrawResponse const& srep, double t0, ComboHit const& chit,  const Tracker* tracker){
        double time_residual_long = TimeResidualLong(doca, srep,  t0,  chit,  tracker);
        double time_residual_trans = TimeResidualTrans(doca);
        return time_residual_trans + time_residual_long;// + hitlen/299 + fltlen/299;
}


