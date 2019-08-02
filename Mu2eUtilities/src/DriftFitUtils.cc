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
TrkPoca DriftFitUtils::GetPOCA(Straw const& straw, double a0, double a1, double b0, double b1) {
		XYZVec track_position(a0,b0,0);
		XYZVec track_direction(a1,b1,1);
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
	
double DriftFitUtils::TimeResidual(Straw const&  straw, double doca, StrawResponse srep){
                StrawId strawid = straw.id();
      		double phi =  0.;
      	        double expectedtime= srep.StrawResponse::driftDistanceToTime(strawid , fabs(doca), phi);
      	        //double time_residual_trans = expectedtime;
      	        double wirehalflen = straw.halfLength();              
      	        double propagation_time = wirehalflen/299.; 
      	        double time_residual_long = expectedtime+propagation_time;
                return time_residual_long;

}

