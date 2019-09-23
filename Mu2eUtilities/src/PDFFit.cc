//Likilhood Functions for minuit fitting

//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>
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
//Tracker Details:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
//Utilities:
#include "Mu2eUtilities/inc/PDFFit.hh"
#include "Mu2eUtilities/inc/DriftFitUtils.hh"
//Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
using namespace mu2e;

#define AVG_VELOCITY 0.065 //mm/ns
#define sigma 10.5 //ns
#define tau 10.7 //ns

float wireradius = 12.5/1000.; //12.5 um in mm 
float strawradius = 2.5; //2.5 mm in mm 

double TimePDFFit::calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1, ComboHit chit)const{
	double doca = DriftFitUtils::GetTestDOCA(straw, a0,a1,b0,b1, chit); 
        return (doca);
}
double TimePDFFit::TimeResidual(Straw straw, double doca, StrawResponse srep, double t0 ,  ComboHit hit)const{
	double tres =  DriftFitUtils::TimeResidual( straw, doca, srep, t0, hit);
	return tres;
}


double TimePDFFit::operator() (const std::vector<double> &x) const
{
  //Name/store parameters in terms of Minuit "x's":
  double a0 = x[0];
  double a1 = x[1];
  double b0 = x[2];
  double b1 = x[3];
  double t0 = x[4]; 
  long double llike = 0;
  
  //Loop through the straws and get DOCA:
  for (size_t i=0;i<this->straws.size();i++){
      double doca = calculate_DOCA(this->straws[i], a0, a1, b0, b1,chits[i]); 
      double time_residual = this->TimeResidual(this->straws[i], doca, this->srep, t0, this->chits[i]);
      double pdf_t =1/sqrt(2*TMath::Pi()*sigma*sigma) * exp(-(time_residual*time_residual)/(2*sigma*sigma));
      //Log Liklihood:
      llike -=log(pdf_t);
      t0 += time_residual/this->straws.size(); 
      
  }

  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0){
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
    }
  }
  
  return llike;
}

int FullFit::Factorial(int k)
{
  if (k == 0){
    return 1;
  }
  int response = 1;
  for (int i=1;i<k;i++){
    response *= i;
    return response;
   }
}

void FullFit::CalculateFullPDF() {
  for (int is=0;is<pdf_sbins;is++){
     double sigma = this->pdf_sigmas[is];
     for (int it0=0;it0<pdf_tbins;it0++){
       //Time residuals?
       double time_gaus = this->pdf_times[it0];
       //Time PDF
       double time_gaussian = 1.0/sqrt(2*TMath::Pi()*sigma*sigma)*exp(-(time_gaus*time_gaus)/(2*sigma*sigma));
      //Use taus:
      for (int itau=0;itau<pdf_taubins;itau++){
        double tau = this->pdf_taus[itau];
        for (int it1=0;it1<pdf_tbins-it0;it1++){
          double time_tau = this->pdf_deltat*it1;
          double val_tau = pow(1/tau,k)*pow(time_tau,k-1)*exp(-time_tau/tau)/(double) Factorial(k-1);
          this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + (it0+it1)] += time_gaussian * val_tau;
        }
      }
    }
  }

  for (int is=0; is < pdf_sbins;is++){
    for (int itau=0;itau<pdf_taubins;itau++){
      double total = 0;
      for (int it=0;it<pdf_tbins;it++){
        total += this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + it];
      }
      //Normalize:
      for (int it=0;it<pdf_tbins;it++){
        this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + it] /= total;
      }
    }
  }
}


