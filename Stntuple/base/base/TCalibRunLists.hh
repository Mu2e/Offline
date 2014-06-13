#ifndef TCALIBRUNLISTS_HH
#define TCALIBRUNLISTS_HH

#include "TObject.h"
#include "TString.h"


class TCalibRunLists: public TObject {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
  Int_t     fCID;
  TString   fCALIB_TABLE;
  TString   fDATA_STATUS;
  Int_t     fCALIB_RUN;
  Int_t     fCALIB_VERSION;
  Int_t     fCREATE_DATE;
  TString   fCREATE_USER;
  Int_t     fTSID;
  Int_t     fUPDATE_DATE;
  TString   fUPDATE_USER;
  TString   fCREATE_HOST;
  TString   fUPDATE_HOST;
  TString   fDAQ_STATUS;
  TString   fALGORITHM;
  Int_t     fPERIOD;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TCalibRunLists();
  ~TCalibRunLists();

  void  Print(const char* Opt = "") const ;

  ClassDef(TCalibRunLists,1)
};
#endif
