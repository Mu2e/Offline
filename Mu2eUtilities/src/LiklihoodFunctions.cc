//Likilhood Functions for minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive parameters from CosmicTrack stored there.

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
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
using namespace mu2e;

namespace LiklihoodFunctions{
        
	EndResult DoFit(CosmicTrackSeed trackseed , StrawResponse srep){
	  
	  std::vector<double> errors(5,0);
	  std::vector<double> seed(5,0);
         
	  EndResult endresult;
	  CosmicTrack cosmictrack = trackseed._track;
          
          //Seed Track using seed fit parameters stored in track info 
	  seed[0] = trackseed._track.FitEquationXYZ.Pos.X();
	  seed[1] = trackseed._track.FitEquationXYZ.Dir.X();
	  seed[2] = trackseed._track.FitEquationXYZ.Pos.Y();
	  seed[3] = trackseed._track.FitEquationXYZ.Dir.Y();//trackseed._track.FitParams.B1;
	  seed[4] = trackseed._t0.t0();
	  /* alignment parameters -  6 DoF:
	  1) Translations:
	  seed[5] = dx;//shifts of straw centre positions relative to truth
	  seed[6] = dy;
          seed[7] = dz;
	  2) Rotations: 
	  seed[8] =  eta;//shift about x
	  seed[9] = zeta; //shift about y
	  seed[10] = epsilon; //shift about z
	  */
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
          TimePDFFit fit(trackseed._straw_chits, trackseed._straws, srep, cosmictrack, constraint_means,constraints,1);
	  //DataFit fit(trackseed._straw_chits, trackseed._straws, srep, cosmictrack, constraint_means,constraints,1);
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
	  int maxfcn = 400;
	  double tolerance = 1000;
          //Define Minimization method as "MIGRAD" (see minuit documentation)
	  ROOT::Minuit2::FunctionMinimum min = migrad(maxfcn, tolerance);
	  ROOT::Minuit2::MnPrint::SetLevel(3);
	  ROOT::Minuit2::operator<<(cout, min);

	  //Will be the results of the fit routine:
	  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
	  double minval = min.Fval();
	  //Define name for parameters
	  endresult.bestfit = results.Params();
	  endresult.bestfiterrors = results.Errors();
          //Store Minuit Covarience if exisits:
	  if(min.HasValidCovariance()) endresult.bestfitcov = min.UserCovariance().Data();
          //Name Parameters:
	  endresult.names.push_back("a0");
	  endresult.names.push_back("a1");
	  endresult.names.push_back("b0");
	  endresult.names.push_back("b1");
	  endresult.names.push_back("t0");
	  
	  //Add best fit results to appropriatly named element:
	  endresult.NLL = minval;
	  for (size_t i=0;i<endresult.names.size();i++){
	    std::cout << i << endresult.names[i] << " : " << endresult.bestfit[i] << " +- " << endresult.bestfiterrors[i] << std::endl;
	    if(endresult.bestfitcov.size() != 0 and i< 4) cout<<"cov "<<endresult.bestfitcov[i]<<endl;
	
	  }
	
	  //Store DOCA results for analysis:
	  if(minval != 0 ){ cosmictrack.minuit_converged = true;} 
	  

	  for(size_t i = 0; i< trackseed._straws.size(); i++){	    
	      double start_doca = fit.calculate_DOCA(trackseed._straws[i],seed[0], seed[1], seed[2], seed[3], trackseed._straw_chits[i]);
	      double start_time_residual = fit.TimeResidual(trackseed._straws[i], start_doca,  srep, seed[4], trackseed._straw_chits[i]);
	      endresult.StartDOCAs.push_back(start_doca);
	      endresult.StartTimeResiduals.push_back(start_time_residual);
	      double end_doca = fit.calculate_DOCA(trackseed._straws[i],endresult.bestfit[0], endresult.bestfit[1], endresult.bestfit[2], endresult.bestfit[3], trackseed._straw_chits[i]);
	      double end_time_residual = fit.TimeResidual(trackseed._straws[i], end_doca,  srep, endresult.bestfit[4], trackseed._straw_chits[i]);
	      endresult.EndTimeResiduals.push_back(end_time_residual);
	      endresult.EndDOCAs.push_back(end_doca);

	      double true_doca = fit.calculate_DOCA(trackseed._straws[i], trackseed._track.TrueFitEquation.Pos.X(), trackseed._track.TrueFitEquation.Dir.X(), trackseed._track.TrueFitEquation.Pos.Y(),trackseed._track.TrueFitEquation.Dir.Y(), trackseed._straw_chits[i]);
	      endresult.TrueDOCAs.push_back(true_doca);
	      
	      
	  }
	
	 return endresult;
 
  }
  
}
