#ifndef ILCPID_H
#define ILCPID_H

/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// particle id probability densities                                        //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

/* $Id: IlcPID.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */


#include <TObject.h>


class IlcPID : public TObject {
 public:
  enum {
    kSPECIES = 5,    // Number of particle species recognized by the PID
    kSPECIESN = 10   // Number of charged+neutral particle species recognized by the DREAM/DREAMTB/PHOS PID
  };
  enum EParticleType {
    kElectron = 0, 
    kMuon = 1, 
    kPion = 2, 
    kKaon = 3, 
    kProton = 4, 
    kPhoton = 5, 
    kPi0 = 6, 
    kNeutron = 7, 
    kKaon0 = 8, 
    kEleCon = 9,
    kUnknown = 10
  };
  static Float_t       ParticleMass(Int_t iType) 
    {return fgkParticleMass[iType];};
  static const char*   ParticleName(Int_t iType) 
    {return fgkParticleName[iType];};
  static Int_t         ParticleCode(Int_t iType) 
    {return fgkParticleCode[iType];};

  IlcPID();
  IlcPID(const Double_t* probDensity, Bool_t charged = kTRUE);
  IlcPID(const Float_t* probDensity, Bool_t charged = kTRUE);
  IlcPID(const IlcPID& pid);
  IlcPID& operator = (const IlcPID& pid);

  Double_t             GetProbability(EParticleType iType,
				      const Double_t* prior) const;
  Double_t             GetProbability(EParticleType iType) const;
  void                 GetProbabilities(Double_t* probabilities,
					const Double_t* prior) const;
  void                 GetProbabilities(Double_t* probabilities) const;
  EParticleType        GetMostProbable(const Double_t* prior) const;
  EParticleType        GetMostProbable() const;
  
  static void          SetPriors(const Double_t* prior,
				 Bool_t charged = kTRUE);
  static void          SetPrior(EParticleType iType, Double_t prior);

  IlcPID&              operator *= (const IlcPID& pid);

 private:

  void                 Init();

  Bool_t               fCharged;                   // flag for charged/neutral
  Double_t             fProbDensity[kSPECIESN];    // probability densities
  static Double_t      fgPrior[kSPECIESN];         // a priori probabilities

  static /*const*/ Float_t fgkParticleMass[kSPECIESN+1]; // particle masses
  static const char*   fgkParticleName[kSPECIESN+1]; // particle names
  static const Int_t   fgkParticleCode[kSPECIESN+1]; // particle codes

  ClassDef(IlcPID, 1)    // particle id probability densities
};


IlcPID operator * (const IlcPID& pid1, const IlcPID& pid2);


#endif
