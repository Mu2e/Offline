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
	 FitResult DoFit(int _diag, CosmicTrackFinderData& trackdata, StrawResponse const& srep , const Tracker* tracker,double max_doca, unsigned int minChits, int MaxLogL, double _gaussTres, double maxTres){
	  
	
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);
          std::vector<double> newseed(5,0);
	  std::vector<double> newerrors(5,0);
	  FitResult FitResult;
	  CosmicTrack cosmictrack = trackdata._tseed._track;
         
          //Seed Gaussian PDF using seed fit parameters stored in track info 
	  seed[0] = trackdata._tseed._track.FitEquationXYZ.Pos.X();
	  seed[1] = trackdata._tseed._track.FitEquationXYZ.Dir.X();
	  seed[2] = trackdata._tseed._track.FitEquationXYZ.Pos.Y();
	  seed[3] = trackdata._tseed._track.FitEquationXYZ.Dir.Y();//trackdata._tseed._track.FitParams.B1;
	  seed[4] = trackdata._tseed._t0.t0();
	  
	  //Seed errors = covarience of parameters in seed fit
	  errors[0] = trackdata._tseed._track.FitParams.Covarience.sigA0; 
	  errors[1] = trackdata._tseed._track.FitParams.Covarience.sigA1;
	  errors[2] = trackdata._tseed._track.FitParams.Covarience.sigB0;
	  errors[3] = trackdata._tseed._track.FitParams.Covarience.sigB1;
	  errors[4] = trackdata._tseed._t0.t0Err();
	  //Constrain to mean = 0 for 4 parameters (T0 might not be so...)
	  std::vector<double> constraint_means(5,0);
	  std::vector<double> constraints(5,0);
	  //Define the PDF used by Minuit:
	 
          GaussianPDFFit fit(trackdata._tseed._straw_chits,  srep, cosmictrack, constraint_means,constraints, _gaussTres, 1,  tracker);
	  
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
          //Store Minuit Covarience if exists:
	  if(min.HasValidCovariance()) FitResult.bestfitcov = min.UserCovariance().Data();

          //Name Parameters:
	  FitResult.names.push_back("a0");
	  FitResult.names.push_back("a1");
	  FitResult.names.push_back("b0");
	  FitResult.names.push_back("b1");
	  FitResult.names.push_back("t0");
	  FitResult.NLL = minval;
	  for (size_t i=0;i<FitResult.names.size();i++){
			if((!isnan(FitResult.NLL) or FitResult.NLL !=0 or FitResult.NLL< MaxLogL) and !isnan(FitResult.bestfit[i])) { 
				cosmictrack.minuit_converged = true;
			}
	  } 
	
	  //Add best fit results to appropriatly named element:
	  
	  if(_diag > 1){
	  for (size_t i=0;i<FitResult.names.size();i++){
	    std::cout << i << FitResult.names[i] << " : " << FitResult.bestfit[i] << " +- " << FitResult.bestfiterrors[i] << std::endl;
	    if(FitResult.bestfitcov.size() != 0) cout<<"cov "<<FitResult.bestfitcov[i]<<endl;
	
	  }
       }
       //Cut on Gaussian results remove "bad" hits
	ComboHitCollection passed_hits;
	
	for(size_t i = 0; i< trackdata._tseed._straw_chits.size(); i++){
	      double gauss_end_doca = fit.calculate_DOCA(trackdata._tseed._straw_chits[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3],  tracker);
	      double gauss_end_time_residual = fit.TimeResidual( gauss_end_doca,  srep, FitResult.bestfit[4], trackdata._tseed._straw_chits[i],  tracker);
	      FitResult.GaussianEndTimeResiduals.push_back(gauss_end_time_residual);
	      FitResult.GaussianEndDOCAs.push_back(gauss_end_doca);
	       if (gauss_end_doca < max_doca and gauss_end_time_residual < maxTres){ 
			passed_hits.push_back(trackdata._tseed._straw_chits[i]);
			
		} else {
			trackdata._tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier); 
			
		}
	}
        
	if(cosmictrack.minuit_converged and passed_hits.size() > minChits) {

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
	 newerrors[4] = trackdata._tseed._t0.t0Err();
	 FullDriftFit fulldriftfit(passed_hits, srep, cosmictrack, constraint_means,constraints,_gaussTres, 1, tracker);
	
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
	 
	  for(size_t i = 0; i< passed_hits.size(); i++){
	     
	      double start_doca = fit.calculate_DOCA(passed_hits[i],seed[0], seed[1], seed[2], seed[3], tracker);
	      double start_time_residual = fit.TimeResidual(start_doca,  srep, seed[4], trackdata._tseed._straw_chits[i], tracker);
	      
		//Store Final Fit DOCA:
	      double end_doca = fit.calculate_DOCA(passed_hits[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3], tracker);
              double ambig = fit.calculate_ambig(passed_hits[i],FitResult.bestfit[0], FitResult.bestfit[1], FitResult.bestfit[2], FitResult.bestfit[3], tracker);
	      double end_time_residual = fit.TimeResidual( end_doca,  srep, FitResult.bestfit[4], passed_hits[i], tracker);
	      
	      FitResult.StartDOCAs.push_back(start_doca);
	      FitResult.StartTimeResiduals.push_back(start_time_residual);
	      FitResult.FullFitEndTimeResiduals.push_back(end_time_residual);
	      FitResult.FullFitEndDOCAs.push_back(end_doca);
	      FitResult.RecoAmbigs.push_back(ambig);
	      
	  }
	  //delete array list to avoid memory leaks:
	  fulldriftfit.DeleteArrays();
     	}
	 return FitResult;
 
  }
  
}
