#ifndef TUSED_SET_HH
#define TUSED_SET_HH

#include "TObject.h"
#include "TString.h"


class TUsedSet: public TObject {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
  Int_t     fJOBSET;
  TString   fPROCESS_NAME;
  Int_t     fPROCESS_RUN;
  TString   fPROC_CALIB_VERSION;
  TString   fDBNAME;
  Int_t     fCREATE_DATE;
  TString   fCREATE_USER;
  Int_t     fTSID;
  Int_t     fUPDATE_DATE;
  TString   fUPDATE_USER;
  Int_t     fPROC_TAG;
  TString   fPROCESS_STATUS;
  Int_t     fUSID;
  TString   fCREATE_HOST;
  TString   fUPDATE_HOST;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TUsedSet();
  ~TUsedSet();

  ClassDef(TUsedSet,1)
};
#endif
