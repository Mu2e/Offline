//Likilhood Functions for minuit fitting
// Author: S. Middleton 
// Date: July 2019


//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Mu2eUtilities/inc/PDFFit.hh"
#include "Mu2eUtilities/inc/DriftUtils.hh"

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

#define AVG_VELOCITY 0.065
float wireradius = 12.5/1000.; //12.5 um in mm 
float strawradius = 2.5; //2.5 mm in mm 


void PDFFit::calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min, double doca_max, bool dist) const
{
  std::cout<<"Get weighted PDF "<<std::endl;
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  std::cout<<"Times size in pdf "<<this->times.size()<<std::endl;
  for (size_t i=0;i<this->times.size();i++){
    if (fabs(this->docas[i]) < doca_min || fabs(this->docas[i]) > doca_max)
      continue;
    std::cout<<"Getting Taus"<<std::endl;
    double hypotenuse = sqrt(pow(this->docas[i],2) + pow(tau * AVG_VELOCITY,2));
    double tau_eff = hypotenuse/AVG_VELOCITY - this->docas[i]/AVG_VELOCITY;
    std::cout<<"Hist Stuff "<<std::endl;
    for (int j=0;j<h->GetNbinsX();j++){
      double time_residual = h->GetXaxis()->GetBinCenter(j+1);
      if (dist){
        time_residual = h->GetXaxis()->GetBinCenter(j+1)/AVG_VELOCITY;
      }

      double pdf_val =1/sqrt(2*TMath::Pi()*x[1]*x[1]) * exp(-(time_residual*time_residual)/(2*x[1]*x[1]))+x[3];
      std::cout<<"Tau "<<tau<<" sigma "<<sigma<<"Background "<< background<<"toffset"<<toffset<<"time residual "<<time_residual<<" PDF "<<pdf_val<<" Tau Eff "<<tau_eff<<std::endl;
    }
  }
}


double PDFFit::operator() (const std::vector<double> &x) const
{
  std::cout<<"In Operator "<<std::endl;
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  double llike = 0;
  
  for (size_t i=0;i<this->docas.size();i++){
    double time_residual = this->times[i]-toffset;//this->TimeResidual(this->docas[i],this->times[i],toffset);
    std::cout<<"t offset "<<toffset<<"time "<<times[i]<<" time res "<<time_residual<<" doca "<<docas[i]<<" taus "<<tau<<" sigma "<<sigma<<std::endl;
    double hypotenuse = sqrt(pow(this->docas[i],2) + pow(tau * AVG_VELOCITY,2));
    std::cout<<"hyptoenuse "<<hypotenuse<<std::endl;
    double tau_eff = hypotenuse/0.0625 - this->docas[i]/0.0625;
    std::cout<<"Tau Eff"<<tau_eff<<std::endl;
    double pdf_val = 1/sqrt(2*TMath::Pi()*x[1]*x[1]) * exp(-(time_residual*time_residual)/(2*x[1]*x[1]))+x[3];
    llike -= log(pdf_val);
  }
  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0){
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
  }
  }
  std::cout<<"Sigma "<<sigma<<" Background "<<background<<"llike"<<llike<<std::endl;
  return llike;
}
