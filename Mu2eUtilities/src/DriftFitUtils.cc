#include "Mu2eUtilities/inc/DriftFitUtils.hh"
#include "Mu2eUtilities/inc/ConvertXYZ.hh"
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
double fltlen = 0;
double hitlen = 0;

void SetFltLen(double length){
	fltlen = length;

}

void SetHitLen(double length){
	hitlen = length;

}


double DriftFitUtils::GetTestDOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit chit) {
	XYZVec track_position(a0,b0,0);
	XYZVec track_direction(a1,b1,1);

	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
	const CLHEP::Hep3Vector& sdir = straw.getDirection();
	
	XYZVec wire_position = ConvertToXYZ(spos);
        XYZVec wire_direction= ConvertToXYZ(sdir);
	
	TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_position,
                track_direction,
                wire_position,
                wire_direction,
                1.e-8);
	
        double dca; 
	dca = PCA.dca();     		
	//ParametricFit::LineToLineDCA(track_position, track_direction, wire_position, wire_direction, dca);
	return dca;
}

TrkPoca DriftFitUtils::GetPOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit chit) {
		//Get the current position and direction vector NB: this is such that the parameter is "z" coordinate and the position is the point in X or Y where z=0 respectively.
		XYZVec track_position(a0,b0,0);
		XYZVec track_direction(a1,b1,1);
		
		track_direction = track_direction.Unit();
	        
                //Get straw details
      		const CLHEP::Hep3Vector& spos = straw.getMidPoint();
      		const CLHEP::Hep3Vector& sdir = straw.getDirection();
		
      		XYZVec wire_position = ConvertToXYZ(spos);
                XYZVec wire_direction= ConvertToXYZ(sdir);
             	
      		HepPoint pointwire = ConvertToHepPoint(wire_position);		
      		HepPoint tpos = ConvertToHepPoint(track_position);
      		Hep3Vector tdir = ConvertToHep3Vector(track_direction);
      		Hep3Vector wdir = ConvertToHep3Vector(wire_direction);
                 

      		double fltlen = (fabs((pointwire.z()-tpos.z())/tdir.z())); 
	        SetFltLen(fltlen); //store
		
	 	double hitlen = fabs(wdir.dot(tpos - pointwire));
	 	SetHitLen(hitlen); //store
		
      		TrkLineTraj wire(pointwire,wdir);
      		TrkLineTraj track(tpos,tdir);
      		
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
     }

double DriftFitUtils::GetAmbig(TrkPoca hitpoca) {
		double newamb = hitpoca.doca() > 0 ? 1 : -1;
      		return newamb;
}
	
double DriftFitUtils::GetPropVelocity(StrawResponse rep, ComboHit chit){
	   	double vprop = 2.0*rep.halfPropV(chit.strawId(),1000.0*chit.energyDep());
	   	return vprop; 
}

double DriftFitUtils::GetPropTime(ComboHit chit, Straw straw, double vprop) {
	    double tprop = 0.;
	    //if(poca.status().success()){
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
                //StrawId strawid = straw.id();
      		//double phi =  hit.phi();
      	        double drift_time= doca/0.065;//srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);
      	        return drift_time;
}
        
double DriftFitUtils::TimeResidualLong(Straw const&  straw, double doca, StrawResponse srep,  double t0, ComboHit hit){
		
                StrawId strawid = straw.id();
      	        double _vprop = 2.0*srep.StrawResponse::halfPropV(strawid,1000.0*hit.energyDep());
      	        double propagation_time = GetPropTime(hit, straw,_vprop);
      	        //double propagation_time = straw.halfLength()/_vprop;
      	        return propagation_time;
}

double DriftFitUtils::TimeResidual(Straw const&  straw, double doca, StrawResponse srep, double t0, ComboHit hit){
		double time_residual_long = TimeResidualLong( straw,  doca, srep,  t0,  hit);
		double time_residual_trans = TimeResidualTrans(straw,doca, srep, t0, hit); 
		return time_residual_trans + time_residual_long + hitlen/299 + fltlen/299;
}

double DriftFitUtils::T0(Straw const&  straw, double doca, StrawResponse srep, double t0, ComboHit hit, double Aver){
		double time_residual_long = TimeResidualLong( straw,  doca, srep,  t0,  hit);
		double time_residual_trans = TimeResidualTrans(straw,doca, srep, t0, hit); 
		Aver += hit.time() - time_residual_trans - time_residual_long;
		return Aver; //hit.time() - time_residual_trans - time_residual_long;
}

