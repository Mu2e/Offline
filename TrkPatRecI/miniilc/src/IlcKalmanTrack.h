#ifndef ILCKALMANTRACK_H
#define ILCKALMANTRACK_H

/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcKalmanTrack.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------------------------
//                          Class IlcKalmanTrack
//      fixed the interface for the derived reconstructed track classes 
//            Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include "IlcExternalTrackParam.h"
#include "IlcLog.h"
#include "IlcPID.h"

class IlcCluster;

class IlcKalmanTrack : public IlcExternalTrackParam {
public:
  IlcKalmanTrack();
  IlcKalmanTrack(const IlcKalmanTrack &t);
  virtual ~IlcKalmanTrack(){};

  void SetLabel(Int_t lab) {fLab=lab;}

  virtual Double_t GetPredictedChi2(const IlcCluster *c) const = 0;
  virtual Bool_t PropagateTo(Double_t xr, Double_t x0, Double_t rho) = 0;
  virtual Bool_t Update(const IlcCluster* c, Double_t chi2, Int_t index) = 0;

  Bool_t   IsSortable() const {return kTRUE;}
  Int_t    GetLabel()   const {return fLab;}
  Double_t GetChi2()    const {return fChi2;}
  Double_t GetMass()    const {return fMass;}
  Int_t    GetNumberOfClusters() const {return fN;}
  virtual Int_t GetClusterIndex(Int_t) const { //reserved for IlcTracker
    IlcWarning("Method must be overloaded !\n");
    return 0;
  } 
  virtual Int_t GetMaxNumberOfLayers() const { //reserved for IlcTracker 
    IlcWarning("Method must be overloaded !\n");
    return 0;
  }

  virtual Double_t GetPIDsignal() const {
    IlcWarning("Method must be overloaded !\n");
    return 0.;
  }

  virtual Int_t Compare(const TObject *) const {return 0;} 

  void GetExternalParameters(Double_t &xr,Double_t p[5]) const {
    xr=GetX();
    for (Int_t i=0; i<5; i++) p[i]=GetParameter()[i];
  }
  void GetExternalCovariance(Double_t cov[15]) const {
    for (Int_t i=0; i<15; i++) cov[i]=GetCovariance()[i];
  }

  static 
  Double_t MeanMaterialBudget(Double_t *start,Double_t *end,Double_t *mparam);

  // Time integration (S.Radomski@gsi.de)
  void StartTimeIntegral();
  void SetIntegratedLength(Double_t l) {fIntegratedLength=l;}
  void SetIntegratedTimes(const Double_t *times);

  Bool_t IsStartedTimeIntegral() const {return fStartTimeIntegral;}
  void AddTimeStep(Double_t length);
  void GetIntegratedTimes(Double_t *times) const;
  Double_t GetIntegratedTime(Int_t pdg) const;
  Double_t GetIntegratedLength() const {return fIntegratedLength;}

  void SetNumberOfClusters(Int_t n) {fN=n;} 

  void SetFakeRatio(Float_t ratio) {fFakeRatio=ratio;}
  Float_t  GetFakeRatio()   const {return fFakeRatio;}
  void SetMass(Double_t mass) {fMass=mass;}

protected:
  void SetChi2(Double_t chi2) {fChi2=chi2;} 

  Int_t fLab;             // track label
  Float_t fFakeRatio;     // fake ratio
  Double_t fChi2;         // total chi2 value for this track
  Double_t fMass;         // mass hypothesis
  Int_t fN;               // number of associated clusters

private:
  // variables for time integration (S.Radomski@gsi.de)
  Bool_t  fStartTimeIntegral;       // indicator wether integrate time
  Double_t fIntegratedTime[IlcPID::kSPECIES];       // integrated time
  Double_t fIntegratedLength;        // integrated length
  
  ClassDef(IlcKalmanTrack,6)    // Reconstructed track
};

#endif


