//Likilhood Functions for minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive parameters from CosmicTrack stored there.
// Author: S. Middleton, based on Tracker Code
// Date: July 2019


//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Mu2eUtilities/inc/LiklihoodFunctions.hh"
#include "Mu2eUtilities/inc/PDFFit.hh"
#include "Mu2eUtilities/inc/ParametricFit.hh"
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
        
	EndResult DoFit(CosmicTrackSeed trackseed , StrawResponse srep){
	  cout<<"================= New End Result==========="<<endl;
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);

	  EndResult endresult;
	  CosmicTrack cosmictrack = trackseed._track;
	  seed[0] = trackseed._track.FitParams.A0;
	  seed[1] = trackseed._track.FitParams.A1;
	  seed[2] = trackseed._track.FitParams.B0;
	  seed[3] = trackseed._track.FitParams.B1;
	  seed[4] = trackseed._t0.t0();
	   
	  errors[0] = trackseed._track.FitParams.Covarience.sigA0;
	  errors[1] = trackseed._track.FitParams.Covarience.sigA1;
	  errors[2] = trackseed._track.FitParams.Covarience.sigB0;
	  errors[3] =trackseed._track.FitParams.Covarience.sigB1;
	  errors[4] = trackseed._t0.t0Err();
	 
	  std::vector<double> constraint_means(5,0);
	  std::vector<double> constraints(5,0);
	  cout<<"Seeds in "<<seed[0]<<" "<<seed[1]<<" "<<seed[2]<<" "<<seed[3]<<endl;
          TimePDFFit fit(trackseed._straw_chits, trackseed._straws, srep, cosmictrack, constraint_means,constraints,1);
        
	  ROOT::Minuit2::MnStrategy mnStrategy(2); 
	  ROOT::Minuit2::MnUserParameters params(seed,errors);
	  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
	  
	  migrad.SetLimits((unsigned) 0, -1500,1500);
	  migrad.SetLimits((unsigned) 1, -1,1);
	  migrad.SetLimits((unsigned) 2, -1500,1500);
	  migrad.SetLimits((unsigned) 3,-1,1);
	  migrad.Fix((unsigned) 4); 
	  int maxfcn = 400;
	  double tolerance = 1.;
	  ROOT::Minuit2::FunctionMinimum min = migrad(maxfcn, tolerance);
	  
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
	  double dif = 0;
	  for (size_t i=0;i<endresult.names.size();i++){
	    dif +=abs(endresult.bestfit[i] - seed[i]);
	    std::cout << i << endresult.names[i] << " : " << endresult.bestfit[i] << " +- " << endresult.bestfiterrors[i] << std::endl;
	  }
	  cout<<"Seeds out "<<seed[0]<<" "<<seed[1]<<" "<<seed[2]<<" "<<seed[3]<<endl;
	  cout<<"diff sum "<<dif<<endl;
	  //if(min.IsValid()){ cosmictrack.minuit_converged = true;}
	  for(size_t i = 0; i< trackseed._straws.size(); i++){
	      //if(cosmictrack.minuit_converged == false) continue;
	      double start_doca = fit.calculate_DOCA(trackseed._straws[i],seed[0], seed[1], seed[2], seed[3], trackseed._straw_chits[i]);
	      double start_time_residual = fit.TimeResidual(trackseed._straws[i], start_doca,  srep, seed[4], trackseed._straw_chits[i]);
	      endresult.StartDOCAs.push_back(start_doca);
	      endresult.StartTimeResiduals.push_back(start_time_residual);
	      double end_doca = fit.calculate_DOCA(trackseed._straws[i],endresult.bestfit[0], endresult.bestfit[1], endresult.bestfit[2], endresult.bestfit[3], trackseed._straw_chits[i]);
	      double end_time_residual = fit.TimeResidual(trackseed._straws[i], end_doca,  srep, endresult.bestfit[4], trackseed._straw_chits[i]);
	      endresult.EndTimeResiduals.push_back(end_time_residual);
	      endresult.EndDOCAs.push_back(end_doca);
	     
	  }
	  
	 return endresult;
 
  }
  
}
