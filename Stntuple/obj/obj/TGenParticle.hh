#ifndef STNTUPLE_TGenParticle
#define STNTUPLE_TGenParticle

#include "TMath.h"
#include "TParticle.h"
#include "TParticlePDG.h"

class GENP_PARTICLE;

class TGenParticle : public TParticle {
					// commonly used PDG codes (antiparticle
					// codes are negative)
  enum {
    kElectronPdgCode =  11,
    kNuePdgCode      =  12,
    kMuonPdgCode     =  13,
    kNumuPdgCode     =  14,
    kTauPdgCode      =  15,
    kNutauPdgCode    =  16,
    kPhotonPdgCode   =  22,
    kPi0PdgCode      = 111,
    kPiPlusPdgCode   = 211,
    kK0sPdgCode      = 310
  };

  float   fProperTime;

public:
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
				// ****** constructors and destructor
  TGenParticle();

  TGenParticle(Int_t ip, Int_t idhep, Int_t istdhep, 
	       Int_t jm1, Int_t jm2, Int_t jd1, Int_t jd2, 
	       Float_t px, Float_t py, Float_t pz, Float_t e,
	       Float_t vx, Float_t vy, Float_t vz, Float_t t,
	       float ProperTime);

  virtual ~TGenParticle();
//-----------------------------------------------------------------------------
// init methods
//-----------------------------------------------------------------------------
  Int_t Init(Int_t   ip , Int_t idhep, Int_t istdhep, 
	     Int_t   jm1, Int_t jm2, Int_t jd1, Int_t jd2, 
	     Float_t px , Float_t py, Float_t pz, Float_t e,
	     Float_t vx , Float_t vy, Float_t vz, Float_t t,
	     float ProperTime);

  Int_t Init(Int_t Ip, const TParticle* P);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t    Number    () const { return GetUniqueID(); }
  Double_t Charge    () const { return ((TParticle*) this)->GetPDG()->Charge(); }
  float    ProperTime() const { return fProperTime; }
//------------------------------------------------------------------------------
//  missing TParticle accessors and setters
//------------------------------------------------------------------------------
  Double_t PolarTheta() { return fPolarTheta; }
  Double_t PolarPhi  () { return fPolarPhi;   }

  void     SetPolarTheta(Double_t Theta) { fPolarTheta = Theta; }
  void     SetPolarPhi  (Double_t Phi  ) { fPolarPhi   = Phi;   }

  Int_t   IsElectron() const { return (TMath::Abs(fPdgCode) == kElectronPdgCode); }
  Int_t   IsMuon    () const { return (TMath::Abs(fPdgCode) == kMuonPdgCode);     }
  Int_t   IsPhoton  () const { return (TMath::Abs(fPdgCode) == kPhotonPdgCode);   }
  Int_t   IsPi0     () const { return (TMath::Abs(fPdgCode) == kPi0PdgCode);      }

  Int_t   IsNeutrino() const { 
    Int_t code = TMath::Abs(fPdgCode);
    if      (code == kNuePdgCode  ) return 1;
    else if (code == kNumuPdgCode ) return 1;
    else if (code == kNutauPdgCode) return 1;
    else                            return 0;
  }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void     Print(Option_t* opt = "") const;

  ClassDef(TGenParticle,1)
};

#endif


