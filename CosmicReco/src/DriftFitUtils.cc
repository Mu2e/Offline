// Author : S Middleton
// Date : 2019
// Purpose : Stores functions associted with Cosmic Drift Fit

#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "TMatrixD.h"

//Tracker Drift Conditions:
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"

//For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

using namespace mu2e;

    double DriftFitUtils::GetTestDOCA(ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker) {
	
	Straw const& straw = tracker->getStraw(chit.strawId());

	XYZVec track_position(a0,b0,0);
	XYZVec track_direction(a1,b1,1);

	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
	const CLHEP::Hep3Vector& sdir = straw.getDirection();

	XYZVec wire_position = Geom::toXYZVec(spos);
	XYZVec wire_direction= Geom::toXYZVec(sdir);

	TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
		track_direction,
		wire_position,
		wire_direction,
		1.e-8);

	double dca; 
	dca = PCA.dca();    
	
	return dca;
}

    double DriftFitUtils::GetTestDOCA(ComboHit const& chit, XYZVec const& track_position, XYZVec const&  track_direction, const Tracker* tracker) {
	
	Straw const& straw = tracker->getStraw(chit.strawId());
	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
	const CLHEP::Hep3Vector& sdir = straw.getDirection();
	
	XYZVec wire_position = Geom::toXYZVec(spos);
        XYZVec wire_direction= Geom::toXYZVec(sdir);
	
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
	XYZVec track_pos(a0,b0,0);
	XYZVec track_dir(a1,b1,1);

 	Hep3Vector td(a1, b1, 1);
	td = td.unit();
	Hep3Vector rperp = td - (td.dot(straw.getDirection()) * straw.getDirection());
	auto const& straw_mp = straw.getMidPoint();
        auto const& wire_dir = straw.getDirection().unit();
	TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_pos, track_dir, Geom::toXYZVec(straw_mp),
		                        Geom::toXYZVec(wire_dir), 1.e-8);

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
	XYZVec track_pos(a0,b0,0);
	XYZVec track_dir(a1,b1,1);

 	Hep3Vector td(a1, b1, 1);
	td = td.unit();
	Hep3Vector rperp = td - (td.dot(straw.getDirection()) * straw.getDirection());
	
	double phi = rperp.theta();
	double drift_distance = _srep.driftTimeToDistance(straw_id, chit.driftTime(),
		                  phi); 
	return drift_distance;
}


    int DriftFitUtils::GetAmbig(ComboHit const& chit, XYZVec const& track_position, XYZVec const&  track_direction, const Tracker* tracker) {
	Straw const& straw = tracker->getStraw(chit.strawId());
	
	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
	const CLHEP::Hep3Vector& sdir = straw.getDirection();
	
	XYZVec wire_position = Geom::toXYZVec(spos);
        XYZVec wire_direction= Geom::toXYZVec(sdir);
	
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
	
	XYZVec track_position(a0,b0,0);
	XYZVec track_direction(a1,b1,1);

	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
	const CLHEP::Hep3Vector& sdir = straw.getDirection();
	
	XYZVec wire_position = Geom::toXYZVec(spos);
        XYZVec wire_direction= Geom::toXYZVec(sdir);
	
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
	switch (chit.driftEnd()) {
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


