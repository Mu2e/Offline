#ifndef TTRIGGER_CONDITIONS_HH
#define TTRIGGER_CONDITIONS_HH

#include "TObject.h"
#include "TString.h"


class TTriggerConditions: public TObject {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  Int_t     fRunNumber;
  TString   fName;
  Int_t     fBit;
  Int_t     fLevel;
  TString   fCreateTime;
  Int_t     fFrontendTime;
  Int_t     fUnprescaled;
  Int_t     fPrescaled;
  Int_t     fLive;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				        // ****** constructors and destructor
  TTriggerConditions();
  ~TTriggerConditions();
					// ****** accessors

  Int_t    RunNumber   () { return fRunNumber; }
  TString& Name        () { return fName; }
  TString& CreateTime  () { return fCreateTime; }
  Int_t    Bit         () { return fBit; }
  Int_t    Level       () { return fLevel; }
  Int_t    Unprescaled () { return fUnprescaled; }
  Int_t    Prescaled   () { return fPrescaled; }
  Int_t    Live        () { return fLive; }
  Int_t    FrontendTime() { return fFrontendTime; }

  ClassDef(TTriggerConditions,1)
};
#endif
