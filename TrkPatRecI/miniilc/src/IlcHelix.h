#ifndef ILCHELIX_H
#define ILCHELIX_H

/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//                          Class IlcHelix
//
//         Origin: Marian Ivanov marian.ivanov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>


class IlcCluster;
class IlcKalmanTrack;
class IlcExternalTrackParam;

class IlcHelix : public TObject {
public:
  IlcHelix();
  IlcHelix(const IlcHelix &t);
  IlcHelix(const IlcKalmanTrack &t);
  IlcHelix(const IlcExternalTrackParam &t);
  IlcHelix(Double_t x[3], Double_t p[3], Double_t charge=1, Double_t conversion=0.);
  virtual ~IlcHelix(){};
  inline void Evaluate(Double_t t, Double_t r[3]);
  void Evaluate(Double_t t, Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
		Double_t gg[3]);     //second derivatives
  void GetMomentum(Double_t phase, Double_t p[4], Double_t conversion=0., Double_t *xr=0);  // return  momentum  
  void  GetAngle(Double_t t1, IlcHelix &h, Double_t t2, Double_t angle[3]);
  inline Double_t GetHelixR(Double_t phase=0);
  inline Double_t GetHelixZ(Double_t phase=0);

  Double_t  GetPhase(Double_t x0, Double_t y0);  //return phase for nearest point  
  Double_t  GetPhaseZ(Double_t z0);               // return phase for given z
  Int_t     GetPhase(Double_t r0, Double_t t[2]);               //return phase for the nearest point
  Int_t     GetRPHIintersections(IlcHelix &h, Double_t phase[2][2], Double_t ri[2], Double_t cut=3.);
  Int_t     GetClosestPhases(IlcHelix &h, Double_t phase[2][2]);
  Double_t  GetPointAngle(IlcHelix &h, Double_t phase[2],const Float_t *vertex);
  Int_t     LinearDCA(IlcHelix &h, Double_t &t1, Double_t &t2, 
		      Double_t &R, Double_t &dist);
  //
  Int_t     ParabolicDCA(IlcHelix&h,  //helixes
			 Double_t &t1, Double_t &t2, 
			 Double_t &R, Double_t &dist, Int_t iter=1);    
  Int_t     ParabolicDCA2(IlcHelix&h,  //helixes
			 Double_t &t1, Double_t &t2, 
			 Double_t &R, Double_t &dist, Double_t err[3], Int_t iter=1);    
  Double_t GetHelix(Int_t i) const{return fHelix[i];}
 public:
  Double_t fHelix[9];    //helix parameters
 private:  
  ClassDef(IlcHelix,1)    // IlcHelix
};

void IlcHelix::Evaluate(Double_t t, Double_t r[3]){
  //
  // calculate poitition at given phase t 
  Double_t phase=fHelix[4]*t+fHelix[2];  
  r[0] = fHelix[5] + TMath::Sin(phase)/fHelix[4];
  r[1] = fHelix[0] - TMath::Cos(phase)/fHelix[4];  
  r[2] = fHelix[1] + fHelix[3]*t;
}

inline Double_t  IlcHelix::GetHelixR(Double_t phase)
{
  Double_t x[3];
  Evaluate(phase,x);
  return TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
}

inline Double_t  IlcHelix::GetHelixZ(Double_t phase)
{
  Double_t x[3];
  Evaluate(phase,x);
  return x[2];
}


#endif


