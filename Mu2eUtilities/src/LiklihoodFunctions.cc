//Likilhood Functions for minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive parameters from CosmicTrack stored there.
// Author: S. Middleton, based on Tracker Code
// Date: July 2019


//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Mu2eUtilities/inc/LiklihoodFunctions.hh"
#include "Mu2eUtilities/inc/PDFFit.hh"
#include "TrackerGeom/inc/Tracker.hh"
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
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>
//Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
using namespace mu2e;


namespace LiklihoodFunctions{
        /*
	double TimeResiduals(double drifttime, double doca, double time, double time_offset){
	        
      	        double time_residual_trans = expectedtime;            
      	        double propagation_time = wirehalflen/299.; 
      	        double time_residual_long = time - expectedtime - propagation_time;
	        return time_residual_long;
	        
	        return 0.;
	}
        *
        double GetPhi(Straw const&  straw, XYZVec track_position, std::vector<XYZVec> TrackAxes, double DOCA){
        	const CLHEP::Hep3Vector& spos = straw.getMidPoint();
        	XYZVec wire_position= ConvertToXYZ(spos);
                XYZVec wire_prime_position(wire_position.Dot(TrackAxes[0]),wire_position.Dot(TrackAxes[1]),wire_position.Dot(TrackAxes[2]));
      		//double phi =  acos(track_direction.Dot(wire_prime_direction))/sqrt((track_direction.Mag2()*wire_prime_direction.Mag2()));
      		double phi = acos(track_position.y() - wire_prime_position.y()/DOCA);
      		std::cout<<"phi"<<pi<<std::endl;
      		return phi;	
        }
        
        
	TrkPoca GetPOCA(Straw const&  straw, std::vector<XYZVec> TrackAxes, XYZVec track_position, XYZVec track_direction) {
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
	 	std::cout<<"flight in poca "<<fltlen<<" flight in poca"<<hitlen<<std::endl;
      		TrkLineTraj wire(pointwire,wdir,-straw.halfLength(),straw.halfLength());
      		TrkLineTraj track(tpos,tdir);
      		TrkPoca hitpoca(track,fltlen,wire,hitlen);
      		return hitpoca;
	}
        */
	EndResult DoFit(CosmicTrackSeed trackseed , StrawResponse srep){
	
	  
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);
	  EndResult endresult;//best, best_errors
	 
	  seed[0] = trackseed._track.FitParams.A0;//10;a0
	  seed[1] = trackseed._track.FitParams.A1;//1;a1
	  seed[2] = trackseed._track.FitParams.B0;//b0
	  seed[3] = trackseed._track.FitParams.B1;//b1
	  seed[4] = trackseed._t0.t0(); //t0
	  errors[0] = trackseed._track.FitParams.Covarience.sigA0;
	  errors[1] = trackseed._track.FitParams.Covarience.sigA1;
	  errors[2] = trackseed._track.FitParams.Covarience.sigB0;
	  errors[3] =trackseed._track.FitParams.Covarience.sigB1;
	  errors[4] = trackseed._t0.t0Err();
	  std::cout <<"Seeds : "<<seed[0]<<"  "<<seed[1]<<"  " <<seed[2]<<" "<<seed[3]<<std::endl;
	  std::cout <<"Errors : "<<errors[0]<<"  "<<errors[1]<<"  " <<errors[2]<<" "<<errors[3]<<std::endl;
	  std::vector<double> constraint_means(5,0);
	  std::vector<double> constraints(5,0);
	 
	 
          TimePDFFit fit(trackseed.hits(), trackseed._straws, srep, constraint_means,constraints,1);
          //PDFFit fit(trackseed.hits(), trackseed._straws, srep, constraint_means,constraints,1);
	  ROOT::Minuit2::MnStrategy mnStrategy(2); 
	  ROOT::Minuit2::MnUserParameters params(seed,errors);
	  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
	  
	  migrad.SetLimits((unsigned) 0, -2000,2000);
	  migrad.SetLimits((unsigned) 1, -1,1);
	  migrad.SetLimits((unsigned) 2, -2000,2000);
	  migrad.SetLimits((unsigned) 3,-1,1);
	  migrad.Fix((unsigned) 4); 
	  //int maxfcn = 10;
	  //double tolerance = 1.;
	  ROOT::Minuit2::FunctionMinimum min = migrad();//maxfcn, tolerance);
	  
	  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
	  
	  double minval = min.Fval();
	
	  endresult.bestfit = results.Params();
	  endresult.bestfiterrors = results.Errors();
	  
	  endresult.names.push_back("a0");
	  endresult.names.push_back("a1");
	  endresult.names.push_back("b0");
	  endresult.names.push_back("b1");
	  endresult.names.push_back("t0");
	  //add best fit results to approprtiate name element:
	  std::cout << "NLL: " << minval << std::endl;
	  cout<<"Is Valid: "<<min.IsValid()<<"N calls "<<min.NFcn()<<endl;
	  for (size_t i=0;i<endresult.names.size();i++){
	    std::cout << i << endresult.names[i] << " : " << endresult.bestfit[i] << " +- " << endresult.bestfiterrors[i] << std::endl;
	  }
	 trackseed._track.MinuitFitParams.A0 = endresult.bestfit[0];
	 trackseed._track.MinuitFitParams.A1 = endresult.bestfit[1];
	 trackseed._track.MinuitFitParams.B0 = endresult.bestfit[2];
	 trackseed._track.MinuitFitParams.B1 = endresult.bestfit[3];
	 return endresult;
 
  }
  
}
