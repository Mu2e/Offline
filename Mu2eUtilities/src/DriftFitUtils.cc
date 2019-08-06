#include "Mu2eUtilities/inc/DriftFitUtils.hh"
#include "Mu2eUtilities/inc/ConvertXYZ.hh"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
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

struct drift_info {
	double drift_time;
	double prop_time;
	double v_prop;

};

TrkPoca DriftFitUtils::GetPOCA(Straw const& straw, double a0, double a1, double b0, double b1) {
		std::cout<<"Initial "<<a0<<b0<<b1<<a1<<std::endl;
		XYZVec track_position(a0,b0,0);
		std::cout<<"Initial "<<track_position<<std::endl;
		XYZVec track_direction(a1,b1,1);
		std::cout<<"Initial "<<track_direction<<std::endl;
		track_direction= track_direction.Unit();
		std::vector<XYZVec> TrackAxes= ParametricFit::GetAxes(track_direction);
		
      		const CLHEP::Hep3Vector& spos = straw.getMidPoint();
      		const CLHEP::Hep3Vector& sdir = straw.getDirection();

      		XYZVec wire_position = ConvertToXYZ(spos);
                XYZVec wire_direction= ConvertToXYZ(sdir);

                XYZVec wire_prime_position(wire_position.Dot(TrackAxes[0]),wire_position.Dot(TrackAxes[1]),wire_position.Dot(TrackAxes[2]));
                XYZVec wire_prime_direction(wire_direction.Dot(TrackAxes[0]), wire_direction.Dot(TrackAxes[1]), wire_direction.Dot(TrackAxes[2]));
     
      		HepPoint pointwire = ConvertToHepPoint(wire_prime_position);
      		HepPoint tpos = ConvertToHepPoint(track_position);
      		Hep3Vector tdir = ConvertToHep3Vector(track_direction);
      		Hep3Vector wdir = ConvertToHep3Vector(wire_prime_direction);
      
      		double fltlen = (fabs((spos.z()-tpos.z())/tdir.z())); //along the track
	 	double hitlen = fabs(wdir.dot(tpos - spos)); //along the straw
	        std::cout<<fltlen<<" "<<hitlen<<std::endl;
      		TrkLineTraj wire(pointwire,wdir,-straw.halfLength(),straw.halfLength());
      		TrkLineTraj track(tpos,tdir);
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
     }
     
TrkPoca DriftFitUtils::GetPOCA(Straw const&  straw, std::vector<XYZVec> TrackAxes, XYZVec track_position, XYZVec track_direction) {
      		const CLHEP::Hep3Vector& spos = straw.getMidPoint();
      		const CLHEP::Hep3Vector& sdir = straw.getDirection();
      		XYZVec wire_position = ConvertToXYZ(spos);
                XYZVec wire_direction= ConvertToXYZ(sdir);
                
                XYZVec wire_prime_position(wire_position.Dot(TrackAxes[0]),wire_position.Dot(TrackAxes[1]),wire_position.Dot(TrackAxes[2]));
                XYZVec wire_prime_direction(wire_direction.Dot(TrackAxes[0]), wire_direction.Dot(TrackAxes[1]), wire_direction.Dot(TrackAxes[2]));
                
      		HepPoint pointwire = ConvertToHepPoint(wire_prime_position);
      		HepPoint tpos = ConvertToHepPoint(track_position);
      		Hep3Vector tdir = ConvertToHep3Vector(track_direction);
      		Hep3Vector wdir = ConvertToHep3Vector(wire_prime_direction);
      		
      		double fltlen = (fabs((spos.z()-tpos.z())/tdir.z())); //along the track
	 	double hitlen = fabs(wdir.dot(tpos - spos)); //along the straw
	 	
      		TrkLineTraj wire(pointwire,wdir,-straw.halfLength(),straw.halfLength());
      		TrkLineTraj track(tpos,tdir);
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
	}
	
double DriftFitUtils::GetDOCA(TrkPoca hitpoca) {
	   
        	double doca = hitpoca.doca();
      		return doca;
	}
/*	
double DriftFitUtils::GetPropVelocity(StrawRespose rep, ComboHit chit){
	 vprop = 2.0*rep.halfPropV(_chit.strawId(),1000.0*chit.energyDep());
	 return vprop; 
}
double DriftFitUtils::GetPropTime(ComboHit chit, Straw straw, TrkPoca poca, double vprop) {

    if( poca().status().success()){
      switch (chit.driftEnd()) {
	case StrawEnd::cal:
	  _stime = (straw().halfLength()+hitLen())/vprop;
	  break;
	case StrawEnd::hv:
	  _stime = (straw().halfLength()-hitLen())/vprop;
	  break;
      }
    } else {
// if we're missing poca information, use time division instead
      switch (_combohit.driftEnd()) {
	case StrawEnd::cal:
	  _stime = (straw().halfLength()+timeDiffDist())/vprop;
	  break;
	case StrawEnd::hv:
	  _stime = (straw().halfLength()-timeDiffDist())/vprop;
	  break;
      }
    }
  }
*/
double DriftFitUtils::TimeResidualTrans(Straw const&  straw, double doca, StrawResponse srep){
                StrawId strawid = straw.id();
      		double phi =  0.;
      	        double drift_time= srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);
      	        return drift_time;
      	        }
        
double DriftFitUtils::TimeResidualLong(Straw const&  straw, double doca, StrawResponse srep){
                StrawId strawid = straw.id();
                double chit_e = 1.;//TODO
      	        double expectedtime = TimeResidualTrans( straw,doca, srep); 
      	        double _vprop = 2.0*srep.StrawResponse::halfPropV(strawid,1000.0*chit_e);
      	        double propagation_time = straw.halfLength()/_vprop;
      	        double time_residual_long = expectedtime+propagation_time;
                return time_residual_long;

}

