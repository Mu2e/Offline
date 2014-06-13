#ifndef TBEAM_CONDITIONS_HH
#define TBEAM_CONDITIONS_HH

#include "TObject.h"
#include "TString.h"
/*---------------------------------------------------------------------
root [2] db.DescribeTable("beam_conditions")
RC_RUNNUMBER                   NUMBER
CREATETIME                     DATE
OUTSIDETEMPERATURE             NUMBER
ACCUMULATORSTACK               NUMBER
TEVSTORE                       NUMBER
TEVENERGY                      NUMBER
TEVCURRENT                     NUMBER
B0LOWBETACURRENT               NUMBER
B0LUMINOSITY                   NUMBER
B0LIVELUMINOSITY               NUMBER
B0INTEGRATEDLUMI               NUMBER
B0INTEGRATEDLIVELUMI           NUMBER
B0PROTONLOSSES                 NUMBER
B0ANTIPROTONLOSSES             NUMBER
SVRAD0                         NUMBER
SVRAD1                         NUMBER
SVRAD2                         NUMBER
SVRAD3                         NUMBER
TEVPROTONS                     NUMBER
TEVANTIPROTONS                 NUMBER
FRONTENDTIME                   NUMBER
-------------------------------------------------------------------------*/

class TBeamConditions: public TObject {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  Int_t     fRunNumber;
  TString   fCreateTime;
  Float_t   fOutsideTemperature;
  Float_t   fAccumulatorStack;
  Int_t     fTevStore;
  Float_t   fTevEnergy;
  Float_t   fTevCurrent;
  Float_t   fB0LowBetaCurrent;
  Float_t   fB0Luminosity;
  Float_t   fB0LiveLuminosity;
  Float_t   fB0IntegratedLumi;
  Float_t   fB0IntegratedLiveLumi;
  Float_t   fB0ProtonLosses;
  Float_t   fB0AntiProtonLosses;
  Float_t   fSvrad0;
  Float_t   fSvrad1;
  Float_t   fSvrad2;
  Float_t   fSvrad3;
  Float_t   fTevProtons;
  Float_t   fTevAntiProtons;
  Int_t     fFrontendTime;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TBeamConditions();
  ~TBeamConditions();
				// ****** accessors
  
  Int_t     RunNumber           () { return fRunNumber; }
  Int_t     TevStore            () { return fTevStore; }
  Float_t   TevEnergy           () { return fTevEnergy; }
  Float_t   B0Luminosity        () { return fB0Luminosity; }
  Float_t   B0LiveLuminosity    () { return fB0LiveLuminosity; }
  Float_t   B0IntegratedLumi    () { return fB0IntegratedLumi; }
  Float_t   B0IntegratedLiveLumi() { return fB0IntegratedLiveLumi; }
  Int_t     FrontendTime        () { return fFrontendTime;     }

  ClassDef(TBeamConditions,1)
};
#endif
