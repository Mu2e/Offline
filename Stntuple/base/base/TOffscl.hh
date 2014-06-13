#ifndef TOFFSCL_HH
#define TOFFSCL_HH

#include "TObject.h"
#include "TString.h"


class TOffscl: public TObject {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
  Int_t     fCID;
  Int_t     fCOMPONENT;
  Float_t   fCALIB;
  TString   fCOMPNAME;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TOffscl();
  ~TOffscl();

  void  Print(const char* Opt = "") const ;

  ClassDef(TOffscl,1)
};
#endif
