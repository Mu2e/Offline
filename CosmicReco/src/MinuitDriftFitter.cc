//Author: S Middleton
//Purpose: calls  minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive parameters from CosmicTrack stored there.

//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "CosmicReco/inc/MinuitDriftFitter.hh"
#include "CosmicReco/inc/PDFFit.hh"
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
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
using namespace mu2e;

namespace MinuitDriftFitter{
        
	FitResult DoFit(int _diag, CosmicTrackSeed trackseed , StrawResponse srep, double max_doca, unsigned int minChits, int MaxLogL, double _gaussTres, double maxTres){
	  
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);
          std::vector<double> newseed(5,0);
	  std::vector<double> newerrors(5,0);
	  FitResult FitResult;
	  CosmicTrack cosmictrack = trackseed._track;
         
          //Seed Gaussian PDF using seed fit parameters stored in track info 
	  seed[0] = trackseed._track.FitEquationXYZ.Pos.X();
	  seed[1] = trackseed._track.FitEquationXYZ.Dir.X();
	  seed[2] = trackseed._track.FitEquationXYZ.Pos.Y();
	  seed[3] = trackseed._track.FitEquationXYZ.Dir.Y();//trackseed._track.FitParams.B1;
	  seed[4] = trackseed._t0.t0();
	  
	  //Seed errors = covarience of parameters in seed fit
	  errors[0] = trackseed._track.FitParams.Covarience.sigA0; 
	  errors[1] = trackseed._track.FitParams.Covarience.sigA1;
	  errors[2] = trackseed._track.FitParams.Covarience.sigB0;
	  errors[3] = trackseed._track.FitParams.Covarience.sigB1;
	  errors[4] = trackseed._t0.t0Err();
	  //Constrain to mean = 0 for 4 parameters (T0 might not be so...)
	  std::vector<double> constraint_means(5,0);
	  std::vector<double> constraints(5,0);
	  //Define the PDF used by Minuit:
	 
          GaussianPDFFit fit(trackseed._straw_chits, trackseed._straws, srep, cosmictrack, constraint_means,constraints, _gaussTres, 1);
	  
          //Initiate Minuit Fit:
	  ROOT::Minuit2::MnStrategy mnStrategy(2); 
	  ROOT::Minuit2::MnUserParameters params(seed,errors);
	  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
	  //Set Limits as tracker dimensions:
	  migrad.SetLimits((signed) 0,-10000, 10000);
	  migrad.SetLimits((signed) 1, -5, 5);
	  migrad.SetLimits((signed) 2,-5000, 5000 ); 
	  migrad.SetLimits((signed) 3, -10,10);
	  migrad.Fix((unsigned) 4); 
	  int maxfcn = MaxLogL;
	  double tolerance = 1000;
          //Define Minimization method as "MIGRAD" (see minuit documentation)
	  ROOT::Minuit2::FunctionMinimum min = migrad(maxfcn, tolerance);
	  if(_diag > 1){
	  	ROOT::Minuit2::MnPrint::SetLevel(3);
	  	ROOT::Minuit2::operator<<(cout, min);
	  }
	  //Will be the results of the fit routine:
	  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
	  double minval = min.Fval();
	  //Define name for parameters
	  FitResult.bestfit = results.Params();
	  FitResult.bestfiterrors = results.Errors();
          //Store Minuit Covarience if exisits:
	  if(min.HasValidCovariance()) FitResult.bestfitcov = min.UserCovariance().Data();

          //Name Parameters:
	  FitResult.names.push_back("a0");
	  FitResult.names.push_back("a1");
	  FitResult.names.push_back("b0");
	  FitResult.names.push_back("b1");
	  FitResult.names.push_back("t0");

	  if(minval != 0 and minval< 100 ){ cosmictrack.minuit_converged = true;} 
	  //Add best fit results to appropriatly named element:
	  FitResult.NLL = minval;
	  if(_diag > 1){
	  for (size_t i=0;i<FitResult.names.size();i++){
	    std::cout << i << FitResult.names[i] << " : " << FitResult.bestfit[i] << " +- " << FitResult.bestfiterrors[i] << std::endl;
	    if(FitResult.bestfitcov.size() != 0 and i< 4) cout<<"cov "<<FitResult.bestfitcov[i]<<endl;
	
	  }
       }
       //Cut on Gaussian results remove "bad" hits
	ComboHitCollection passed_hits;
	std::vector<Straw> passed_straws;
	for(size_t i = 0; i< trackseed._straws.size(); i++){
	      double gauss_end_doca = fit.calculate_DOCA(trackseed._straws[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3]);
	      double gauss_end_time_residual = fit.TimeResidual(trackseed._straws[i], gauss_end_doca,  srep, FitResult.bestfit[4], trackseed._straw_chits[i]);
	      FitResult.GaussianEndTimeResiduals.push_back(gauss_end_time_residual);
	      FitResult.GaussianEndDOCAs.push_back(gauss_end_doca);
	       if (gauss_end_doca < max_doca and gauss_end_time_residual < maxTres){ 
			passed_hits.push_back(trackseed._straw_chits[i]);
			passed_hits.push_back(trackseed._straw_chits[i]);
			
		} 
	}
       
	if(cosmictrack.minuit_converged ==true and trackseed._straw_chits.size() > minChits) {

	 //Now run Full Fit:
	 newseed[0] = FitResult.bestfit[0];
	 newseed[1] = FitResult.bestfit[1];
	 newseed[2] = FitResult.bestfit[2];
	 newseed[3] = FitResult.bestfit[3];
	 newseed[4] = FitResult.bestfit[4];

	 newerrors[0] = FitResult.bestfiterrors[0]; 
	 newerrors[1] = FitResult.bestfiterrors[1];
	 newerrors[2] = FitResult.bestfiterrors[2];
	 newerrors[3] = FitResult.bestfiterrors[3];
	 newerrors[4] = trackseed._t0.t0Err();
	 FullDriftFit fulldriftfit(passed_hits, passed_straws, srep, cosmictrack, constraint_means,constraints,_gaussTres, 1);
	
	 ROOT::Minuit2::MnUserParameters newparams(newseed,newerrors);
	 ROOT::Minuit2::MnMigrad newmigrad(fulldriftfit, newparams,mnStrategy);
	 //Set Limits as tracker dimensions:
	 newmigrad.SetLimits((signed) 0,-10000, 10000);
	 newmigrad.SetLimits((signed) 1, -5, 5);
	 newmigrad.SetLimits((signed) 2,-5000, 5000 ); 
	 newmigrad.SetLimits((signed) 3, -10,10);
	 newmigrad.Fix((unsigned) 4); 
	 
         //Define Minimization method as "MIGRAD" (see minuit documentation)
	 min = newmigrad(MaxLogL, tolerance);
	
	 //Will be the results of the fit routine:
	 results = min.UserParameters();
	 minval = min.Fval();
	 //Define name for parameters
	 FitResult.bestfit = results.Params();
	 FitResult.bestfiterrors = results.Errors();
	 
	  for(size_t i = 0; i< trackseed._straws.size(); i++){
		//Store Init DOCA	    
	      double start_doca = fit.calculate_DOCA(trackseed._straws[i],seed[0], seed[1], seed[2], seed[3]);
	      double start_time_residual = fit.TimeResidual(trackseed._straws[i], start_doca,  srep, seed[4], trackseed._straw_chits[i]);
	      
		//Store Final Fit DOCA
	      double end_doca = fit.calculate_DOCA(trackseed._straws[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3]);
              double ambig = fit.calculate_ambig(trackseed._straws[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3]);
	      double end_time_residual = fit.TimeResidual(trackseed._straws[i], end_doca,  srep, FitResult.bestfit[4], trackseed._straw_chits[i]);
	      
	      FitResult.StartDOCAs.push_back(start_doca);
	      FitResult.StartTimeResiduals.push_back(start_time_residual);
	      FitResult.FullFitEndTimeResiduals.push_back(end_time_residual);
	      FitResult.FullFitEndDOCAs.push_back(end_doca);
	      FitResult.RecoAmbigs.push_back(ambig);
	     
	  }
	  fulldriftfit.DeleteArrays();
     	}
	 return FitResult;
 
  }
  
}
