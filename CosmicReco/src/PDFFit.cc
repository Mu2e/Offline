//PDF Functions for minuit fitting

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
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
//Tracker Details:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
//Utilities:
#include "CosmicReco/inc/PDFFit.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"
//Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
using namespace mu2e;

//For Gaussian:
const double drift_v = 0.065; //mm/ns TODO fix this in the code
const double sigma = 1.0; //ns (for full fit only) , if gaussian the sigma is set to 24ns in the FCL parametets
const double tau = 10.7; //ns
//For Full fit:
const int N_tbins = 500;
const int N_taubins = 50;
const int N_sbins = 50;

float wireradius = 12.5/1000.; //12.5 um in mm 
float strawradius = 2.5; //2.5 mm in mm 

FullDriftFit::FullDriftFit(ComboHitCollection _chits, StrawResponse const& _srep , CosmicTrack _track, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double _sigma_t, int _k, const Tracker* _tracker) : GaussianPDFFit(_chits,  _srep, _track,  _constraint_means, _constraints, _sigma_t, _k, _tracker)
{
  //create pdf bins using pre defined numbers:
  pdf_times = new double[N_tbins];
  pdf_taus = new double[N_taubins];
  pdf_sigmas = new double[N_sbins];
  pdf = new double[N_tbins*N_taubins*N_sbins];
  //set min and maxes:

  Min_t = 0.1;
  Min_tau = 5;
  Min_s = -1;//10;(1 for upto  sig = 5)
  Max_t = 1000;
  Max_tau = 15.;
  Max_s = 30.;
  k=1;
  //bin widths
  delta_T = (Max_t-Min_t)/((double) N_tbins);
  delta_S = (Max_s-Min_s)/((double) N_sbins);
  delta_Tau = (Max_tau-Min_tau)/((double) N_taubins);
  //Fill:
  for (int i=0;i<N_tbins;i++){
    pdf_times[i] = Min_t + delta_T * i;
  }
  for (int i=0;i<N_sbins;i++){
    pdf_sigmas[i] = Min_s + delta_S * i;
  }
  for (int i=0;i<N_taubins;i++){
    pdf_taus[i] = Min_tau + delta_Tau * i;
  }
  this->CalculateFullPDF();
};
/*-------Factorial ------//
calculates factorial for deominator
//-----------------------*/
int FullDriftFit::Factorial(int k)
{
  if (k == 0)
    return 1;
  int response = 1;
  for (int i=1;i<k;i++)

    response *= i;
  return response;
}
/* ------- Calc. Full PDF--------/
Fills PDF bins
//----------------------------*/
void FullDriftFit::CalculateFullPDF() {
  
  for (int is=0;is<N_sbins;is++){
     double sigma = this->pdf_sigmas[is];
     for (int it0=0;it0<N_tbins;it0++){
       
       double time_gaus = this->pdf_times[it0];
       double time_gaussian = 1.0/sqrt(2*TMath::Pi()*sigma*sigma)*exp(-(time_gaus*time_gaus)/(2*sigma*sigma));
     
      for (int itau=0;itau<N_taubins;itau++){
        double tau = this->pdf_taus[itau];
        for (int it1=0;it1<N_tbins-it0;it1++){
          double time_tau = this->delta_T*it1;
          
          double val_tau = pow(1/tau,k)*pow(time_tau,k-1)*exp(-time_tau/tau)/(double) Factorial(k-1);
          this->pdf[is * N_taubins * N_tbins + itau * N_tbins + (it0+it1)] += time_gaussian * val_tau;
         
        }
      }
    }
  }

  for (int is=0; is < N_sbins;is++){
    for (int itau=0;itau<N_taubins;itau++){
      double total = 0;
      for (int it=0;it<N_tbins;it++){
        total += this->pdf[is * N_taubins * N_tbins + itau * N_tbins + it];
      }
      
      for (int it=0;it<N_tbins;it++){
        this->pdf[(is*N_taubins*N_tbins)+(itau*N_tbins) + it] /= total;
       
      }
    }
  }
}

void FullDriftFit::DeleteArrays() const{
    
    delete []  pdf_sigmas;
    delete []  pdf_taus;
    delete []  pdf_times;
    delete []  pdf;

}

double FullDriftFit::InterpolatePDF(double time_residual, double sigma, double tau) const
{

  int bin_s = (sigma - this->Min_s)/(this->delta_S);
  if (bin_s >= N_sbins-1)
    bin_s = N_sbins-2;
  if (bin_s < 0)
    bin_s = 0;
 
  double s_d = (sigma - (bin_s*this->delta_S+this->Min_s))/(this->delta_S);
  int bin_tau = (tau - this->Min_tau)/(this->delta_Tau);

  if (bin_tau >= N_taubins-1)
    bin_tau = N_taubins-2;
  double tau_d = (tau - (bin_tau*this->delta_Tau+this->Min_tau))/(this->delta_Tau);

  int bin_t = (time_residual - this->Min_t)/(this->delta_T);
 
  if (bin_t >= N_tbins-1)
    bin_t = N_tbins-2;
  if (bin_t < 0)
    bin_t = 0;
  double t_d = (time_residual - (bin_t*this->delta_T+this->Min_t))/(this->delta_T);
 
  double pdf_val000 = this->pdf[(bin_s+0) * N_taubins * N_tbins + (bin_tau+0) * N_tbins + (bin_t+0)];

  double pdf_val100 = this->pdf[(bin_s+1) * N_taubins * N_tbins + (bin_tau+0) * N_tbins + (bin_t+0)];

  double pdf_val001 = this->pdf[(bin_s+0) * N_taubins * N_tbins + (bin_tau+0) * N_tbins + (bin_t+1)];
 
  double pdf_val101 = this->pdf[(bin_s+1) * N_taubins * N_tbins + (bin_tau+0) * N_tbins + (bin_t+1)];
  
  double pdf_val010 = this->pdf[(bin_s+0) * N_taubins * N_tbins + (bin_tau+1) * N_tbins + (bin_t+0)];
  
  double pdf_val110 = this->pdf[(bin_s+1) * N_taubins * N_tbins + (bin_tau+1) * N_tbins + (bin_t+0)];
  
  double pdf_val011 = this->pdf[(bin_s+0) * N_taubins * N_tbins + (bin_tau+1) * N_tbins + (bin_t+1)];
 
  double pdf_val111 = this->pdf[(bin_s+1) * N_taubins * N_tbins + (bin_tau+1) * N_tbins + (bin_t+1)];
  
  double pdf_val00 = pdf_val000*(1-s_d) + pdf_val100*s_d;
  double pdf_val01 = pdf_val001*(1-s_d) + pdf_val101*s_d;
  double pdf_val10 = pdf_val010*(1-s_d) + pdf_val110*s_d;
  double pdf_val11 = pdf_val011*(1-s_d) + pdf_val111*s_d;
  double pdf_val0 = pdf_val00*(1-tau_d) + pdf_val10*tau_d;
  double pdf_val1 = pdf_val01*(1-tau_d) + pdf_val11*tau_d;
  double pdf_val = pdf_val0*(1-t_d) + pdf_val1*t_d;
  
  return pdf_val;
 
}

// This 3 functions talk to the drift util:
double GaussianPDFFit::calculate_DOCA(ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker)const{
	double doca = DriftFitUtils::GetTestDOCA(chit, a0,a1,b0,b1, tracker); 
        return (doca);
}

double GaussianPDFFit::calculate_ambig(ComboHit const& chit, double a0, double a1, double b0, double b1, const Tracker* tracker)const{
	double ambig = DriftFitUtils::GetAmbig(chit, a0,a1,b0,b1, tracker); 
        return (ambig);
}

double GaussianPDFFit::TimeResidual(double doca, StrawResponse const& srep , double t0 ,  ComboHit const& hit, const Tracker* tracker)const{
	double tres =  DriftFitUtils::TimeResidual(doca, srep, t0, hit, tracker);
	return (tres);
}

//Gaussian PDF function:
double GaussianPDFFit::operator() (const std::vector<double> &x) const
{
  //Name/store parameters in terms of Minuit "x's":
  double a0 = x[0];
  double a1 = x[1];
  double b0 = x[2];
  double b1 = x[3];
  double t0 = x[4]; 
  long double llike = 0;
  
  //Loop through the straws and get DOCA:
  for (size_t i=0;i<this->chits.size();i++){
      double doca = calculate_DOCA(this->chits[i], a0, a1, b0, b1,  this->tracker); 
      double time_residual = this->TimeResidual(doca, this->srep, t0,this->chits[i], this->tracker);
      double pdf_t = 1/sqrt(2*TMath::Pi()*this->sigma_t*this->sigma_t) * exp(-(time_residual*time_residual)/(2*this->sigma_t*this->sigma_t));
      //Log Liklihood:
      llike -=log(pdf_t);
      t0 += time_residual/this->chits.size(); 
      
  }

  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0){
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
    }
  }
  
  return llike;
}

//Full fit PDF function:
double FullDriftFit::operator() (const std::vector<double> &x) const
{
  
  if (sigma > Max_s)
    return 1e10;
  double a0 = x[0];
  double a1 = x[1];
  double b0 = x[2];
  double b1 = x[3];
  double t0 = x[4]; 
  long double llike = 0;
  
  for (size_t i=0;i<this->chits.size();i++){
    double doca = calculate_DOCA(this->chits[i], a0, a1, b0, b1, this->tracker); 
    double doca_penalty = 1; 
    double time_residual = this->TimeResidual(doca, this->srep, t0, this->chits[i], this->tracker);
    
    if(time_residual>40){
	time_residual = 40;
    }
    
    double hypotenuse = sqrt(pow(doca,2) + pow((tau * 0.0625),2));
   
    double tau_eff = (hypotenuse/0.0625) - (doca/0.0625);
   
    double pdf_val = this->InterpolatePDF(time_residual,sigma,tau_eff);
    pdf_val *= doca_penalty;//unused

    if (pdf_val < 1e-3){
      pdf_val = 1e-3;
    }
    llike -= log(pdf_val);
    
  }

  for (int i=0;i<7;i++){
    if (this->constraints[i] > 0)
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
  }

  return llike;
}

