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
	
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);
	  std::vector<double> ds;
	  EndResult endresult;//best, best_errors
	  CosmicTrack cosmictrack = trackseed._track;
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
	 
	  std::vector<double> constraint_means(5,0);
	  std::vector<double> constraints(5,0);
	 
          TimePDFFit fit(trackseed.hits(), trackseed._straws, srep, cosmictrack, ds, constraint_means,constraints,1);
          //PDFFit fit(trackseed.hits(), trackseed._straws, srep, constraint_means,constraints,1);
	  ROOT::Minuit2::MnStrategy mnStrategy(2); 
	  ROOT::Minuit2::MnUserParameters params(seed,errors);
	  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
	  
	  migrad.SetLimits((unsigned) 0, -1500,1500);
	  migrad.SetLimits((unsigned) 1, -1,1);
	  migrad.SetLimits((unsigned) 2, -1500,1500);
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
	  std::cout <<"Seeds : "<<seed[0]<<"  "<<seed[1]<<"  " <<seed[2]<<" "<<seed[3]<<std::endl;
	  if(min.IsValid()){ cosmictrack.minuit_converged = true;}
	  //cout<<"size of docas "<<fit.docas.size()<<endl;
	  
	  //cout<<"start doca "<<fit.docas[0]<<endl;
	  //cout<<"storing doca "<<fit.docas[fit.docas.size()-1]<<endl;
	  //cosmictrack.DriftDiag.StartDOCAs = fit.docas;
	  cosmictrack.DriftDiag.EndDOCAs = ds;//fit.GetDOCAList();
	  //cout<<"docas "<<fit.GetDOCAList()[0]<<"...."<<endl;
	  trackseed._track = cosmictrack;
	  
	 return endresult;
 
  }
  
}
