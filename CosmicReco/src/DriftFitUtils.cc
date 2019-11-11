#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "TMatrixD.h"

//Tracker Drift Conditions:
#include "TrackerGeom/inc/Tracker.hh"
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

double DriftFitUtils::GetTestDOCA(Straw const& straw, double a0, double a1, double b0, double b1) {
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

int DriftFitUtils::GetAmbig(Straw const& straw, double a0, double a1, double b0, double b1) {
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
	
double DriftFitUtils::GetPropVelocity(StrawResponse rep, ComboHit chit){
	   	double vprop = 2.0*rep.halfPropV(chit.strawId(),1000.0*chit.energyDep());
	   	return vprop; 
}

double DriftFitUtils::GetPropTime(ComboHit chit, Straw straw, double vprop) {
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


double DriftFitUtils::TimeResidualTrans(Straw const&  straw, double doca, StrawResponse srep, double t0, ComboHit hit){ 
                
      	        double drift_time= doca/0.065;//srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);//TODO fix the drift vel
      	        return drift_time;
}
        
double DriftFitUtils::TimeResidualLong(Straw const&  straw, double doca, StrawResponse srep,  double t0, ComboHit hit){
		
                StrawId strawid = straw.id();
      	        double _vprop = 2.0*srep.StrawResponse::halfPropV(strawid,1000.0*hit.energyDep());
      	        double propagation_time = GetPropTime(hit, straw,_vprop);
      	        return propagation_time;
}

double DriftFitUtils::TimeResidual(Straw const&  straw, double doca, StrawResponse srep, double t0, ComboHit hit){
		double time_residual_long = TimeResidualLong( straw,  doca, srep,  t0,  hit);
		double time_residual_trans = TimeResidualTrans(straw,doca, srep, t0, hit); 
		return time_residual_trans + time_residual_long;// + hitlen/299 + fltlen/299;
}


