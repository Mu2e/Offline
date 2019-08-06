//Likilhood Functions for minuit fitting
// Author: S. Middleton 
// Date: July 2019
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

#define AVG_VELOCITY 0.065
float wireradius = 12.5/1000.; //12.5 um in mm 
float strawradius = 2.5; //2.5 mm in mm 

double TimePDFFit::calculate_DOCA(Straw const& straw, double a0, double a1, double b0, double b1)const{
	TrkPoca poca = DriftFitUtils::GetPOCA(straw, a0,a1,b0,b1);
	return poca.doca();
}
double TimePDFFit::TimeResidual(Straw straw, double doca, double time, StrawResponse srep)const{
	std::cout<<"time "<<time<<std::endl;
	double tres = time - doca/0.0625;
	DriftFitUtils::TimeResidualLong( straw, doca, srep);
	return tres;
}
double TimePDFFit::operator() (const std::vector<double> &x) const
{

  double a0 = x[0];
  double a1 = x[1];
  double b0 = x[2];
  double b1 = x[3];
  double t0 = x[4];
  //tau = 10, 
  double sigma = 1.;
  double llike = 0;
  std::cout<<"times size "<<this->combohit_times.size()<<std::endl;
  for (size_t i=0;i<this->combohit_times.size();i++){
      
      double doca = calculate_DOCA(this->straws[i], x[0],x[1], x[2], x[3]);
      std::cout<<"DOCA"<<doca<<std::endl;
      double time_residual = this->TimeResidual(this->straws[i], doca,this->combohit_times[i],this->srep);
      double pdf_t =1/sqrt(2*TMath::Pi()*sigma*sigma) * exp(-(time_residual*time_residual)/(2*sigma*sigma));
      llike -=log(pdf_t);
      
  }
  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0){
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
    }
  }
  std::cout<<"a0 "<<a0<<" a1 "<<a1<<" b0 "<<b0<<" b1 "<<b1<<" t0 "<<t0<<std::endl;
  return llike;
}

double PDFFit::operator() (const std::vector<double> &x) const
{
  
  double a0 = x[0];
  double a1 = x[1];
  double b0 = x[2];
  double b1 = x[3];
  double t0 = x[4];
  //tau = 10, sigma = 1
  double llike = 0;
  //for (size_t i=0;i<this->docas.size();i++){
  for (size_t i=0;i<this->combohit_times.size();i++){
      double pdf_x = 1/sqrt(2*TMath::Pi()*this->errorsX[i]*this->errorsX[i]) * exp(-((this->combohit_positions[i].x() - x[0]+x[1]*combohit_positions[i].x())*(this->combohit_positions[i].x() - x[0]+x[1]*combohit_positions[i].x()))/(this->errorsX[i]*this->errorsX[i]));
      double pdf_y = 1/sqrt(2*TMath::Pi()*this->errorsY[i]*this->errorsY[i]) * exp(-((this->combohit_positions[i].y() - x[2]+x[3]*combohit_positions[i].y())*(this->combohit_positions[i].y() - x[2]+x[3]*combohit_positions[i].y()))/(this->errorsY[i]*this->errorsY[i]));
     
      double pdf_val = pdf_x*pdf_y;
      llike -=log(pdf_val);
      
  }
  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0){
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
    }
  }
  std::cout<<"a0 "<<a0<<" a1 "<<a1<<" b0 "<<b0<<" b1 "<<b1<<" t0 "<<t0<<std::endl;
 
  return llike;
}
