#include "Mu2eUtilities/inc/DriftFitUtils.hh"
#include "Mu2eUtilities/inc/ConvertXYZ.hh"
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

/*
The old formalism (NB WILL DELETE ONCE WORKING)
TrkPoca DriftFitUtils::GetPOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit chit) {
		//Get the current position and direction vectora
		XYZVec track_position(a0,b0,0);
		XYZVec track_direction(a1,b1,1);
		track_direction= track_direction.Unit();
                //Get the current co-ordinate system:
		TrackAxes MinuitCoordSystem = ParametricFit::GetTrackAxes(track_direction);
		
                //Get straw details
      		const CLHEP::Hep3Vector& spos = straw.getMidPoint();
      		const CLHEP::Hep3Vector& sdir = straw.getDirection();
      		//Get wire details:
 		//CLHEP::Hep3Vector wpos = //spos + chit.wireDist()*sdir;//TODO - should be just spos
      		XYZVec wire_position = ConvertToXYZ(spos);
                XYZVec wire_direction= ConvertToXYZ(sdir);
             
                XYZVec wire_prime_position(wire_position.Dot(MinuitCoordSystem._XDoublePrime),wire_position.Dot(MinuitCoordSystem._YDoublePrime),wire_position.Dot(MinuitCoordSystem._ZPrime));
                XYZVec wire_prime_direction(wire_direction.Dot(MinuitCoordSystem._XDoublePrime), wire_direction.Dot(MinuitCoordSystem._YDoublePrime), wire_direction.Dot(MinuitCoordSystem._ZPrime));
      	     
      		HepPoint pointwire = ConvertToHepPoint(wire_prime_position);
      		HepPoint tpos = ConvertToHepPoint(track_position);
      		Hep3Vector tdir = ConvertToHep3Vector(track_direction);
      		Hep3Vector wdir = ConvertToHep3Vector(wire_prime_direction);
                
      		double fltlen = (fabs((spos.z()-tpos.z())/tdir.z())); //along the track = spos.z() in  prime frame
	        SetFltLen(fltlen);
	 	double hitlen = fabs(wdir.dot(tpos - spos)); //distance from track to wire centre, i.e dist along the straw
	 	SetHitLen(hitlen);
      		TrkLineTraj wire(pointwire,wdir);//NOTE: Print statmenet comes from here, chit.wireDist()-chit.wireRes(),chit.wireDist()+chit.wireRes());
      		TrkLineTraj track(tpos,tdir);
      		
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
     }


*/
TrkPoca DriftFitUtils::GetPOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit chit) {
		//Get the current position and direction vector NB: this is such that the parameter is "z" coordinate and the position is the point in X or Y where z=0 respectively.
		XYZVec track_position(a0,b0,0);
		XYZVec track_direction(a1,b1,1);
		cout<<"Debugging "<<track_position<<" "<<track_direction<<endl;
		track_direction = track_direction.Unit();
	        cout<<" unit "<<track_direction<<endl;
                
                //Get straw details
      		const CLHEP::Hep3Vector& spos = straw.getMidPoint();
      		const CLHEP::Hep3Vector& sdir = straw.getDirection();
		cout<<" straws "<<spos<<" "<<sdir<<endl;
      		
      		XYZVec wire_position = ConvertToXYZ(spos);
                XYZVec wire_direction= ConvertToXYZ(sdir);
             	cout<<" wire "<<wire_position<<" "<<wire_direction<<endl;
	        cout<<" dir unit "<<wire_direction<<endl;

      		HepPoint pointwire = ConvertToHepPoint(wire_position);		
      		HepPoint tpos = ConvertToHepPoint(track_position);
      		Hep3Vector tdir = ConvertToHep3Vector(track_direction);
      		Hep3Vector wdir = ConvertToHep3Vector(wire_direction);
                 

      		double fltlen = (fabs((pointwire.z()-tpos.z())/tdir.z())); 
	        SetFltLen(fltlen); //store
		cout<<"FltLen "<<fltlen<<endl;
	 	double hitlen = fabs(wdir.dot(tpos - pointwire));
	 	SetHitLen(hitlen); //store
		cout<<"HitLen "<<hitlen<<endl;

      		TrkLineTraj wire(pointwire,wdir);
      		TrkLineTraj track(tpos,tdir);
      		
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
     }

double DriftFitUtils::GetDOCA(TrkPoca hitpoca) {
        	double doca = hitpoca.doca();
      		return doca;
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
	        
                StrawId strawid = straw.id();
      		double phi =  hit.phi();
      	        double drift_time= srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);
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
		//std::cout<<hit.time()<<"trans "<<time_residual_trans<<" long time "<<time_residual_long<<std::endl;
		return time_residual_trans + time_residual_long + hitlen/299 + fltlen/299;
}

double DriftFitUtils::T0(Straw const&  straw, double doca, StrawResponse srep, double t0, ComboHit hit, double Aver){
		double time_residual_long = TimeResidualLong( straw,  doca, srep,  t0,  hit);
		double time_residual_trans = TimeResidualTrans(straw,doca, srep, t0, hit); 
		Aver += hit.time() - time_residual_trans - time_residual_long;
		return Aver; //hit.time() - time_residual_trans - time_residual_long;
}

