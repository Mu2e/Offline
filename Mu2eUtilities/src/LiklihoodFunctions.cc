//Likilhood Functions for minuit fitting
// Author: S. Middleton, based on Tracker Code
// Date: July 2019


//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Mu2eUtilities/inc/LiklihoodFunctions.hh"
#include "Mu2eUtilities/inc/DriftUtils.hh"
#include "Mu2eUtilities/inc/PDFFit.hh"

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

	EndResult DoFit(std::vector<double> times, std::vector<double> docas, std::vector<double> time_residuals){
	  
	  EndResult endresults;//best, best_errors
	  double voltage = 1425.;
	  std::vector<double> seed(4,0);
	  std::vector<double> errors(4,0);
	  seed[0] = 10;
	  seed[1] = 1;
	  seed[2] = 0;
	  seed[3] = 0;
	  errors[0] = 0.1;
	  errors[1] = 0.1;
	  errors[2] = 0.1;
	  errors[3] = 0.1;

	  std::vector<double> constraint_means(4,0);
	  std::vector<double> constraints(4,0);
	  std::cout << times.size() << " total events" << std::endl;
	  PDFFit fit(docas,times,time_residuals, constraint_means,constraints,1,voltage);
          
	  ROOT::Minuit2::MnStrategy mnStrategy(2); 
	  ROOT::Minuit2::MnUserParameters params(seed,errors);
	  std::cout<<"Starting Minuit "<<std::endl;
	  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);
	  std::cout<<"Setting Limits "<<std::endl;//caused SEG FAULT:
	  //migrad.SetLimits((unsigned) 0, fit.pdf_mintau,fit.pdf_maxtau);
	  //migrad.SetLimits((unsigned) 1, fit.pdf_mins,fit.pdf_maxs);
	  
	  migrad.Fix((unsigned) 3);
          std::cout<<"Function min "<<std::endl;
	  ROOT::Minuit2::FunctionMinimum min = migrad();
	  std::cout<<"User params "<<std::endl;
	  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
	  std::cout<<"Getting minval "<<std::endl;
	  double minval = min.Fval();
	  //define best results as the outcome params from minuit
	  std::cout<<"Best Fit "<<std::endl;
	  endresult.bestfit = results.Params();
	  endresult.bestfiterrors = results.Errors();
	  std::cout<<" Printing Out "<<std::endl;
	 
	  endresult.names.push_back("tau (ns)");
	  endresult.names.push_back("sigma (ns)");
	  endresult.names.push_back("time offset (ns)");
	  endresult.names.push_back("background (frac)");
	  //add best fit results to approprtiate name element:
	  std::cout << "NLL: " << minval << std::endl;
	  for (size_t i=0;i<names.size();i++){
	    std::cout << i << endresult.names[i] << " : " << endresult.bestfit[i] << " +- " << endresult.bestfiterrors[i] << std::endl;
	  }
	 return endresult;
 
  }
  
}
