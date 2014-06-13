#ifndef TSET_RUN_MAPS_HH
#define TSET_RUN_MAPS_HH

#include "TObject.h"
#include "TString.h"


class TSetRunMaps: public TObject {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
  TString   fPROCESS_NAME;
  Int_t     fCID;
  Int_t     fJOBSET;
  Int_t     fCREATE_DATE;
  TString   fCREATE_USER;
  Int_t     fTSID;
  Int_t     fUPDATE_DATE;
  TString   fUPDATE_USER;
  TString   fCREATE_HOST;
  TString   fUPDATE_HOST;
  Int_t     fFINISHED;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TSetRunMaps();
  ~TSetRunMaps();

  void  Print(const char* Opt = "") const ;

  ClassDef(TSetRunMaps,1)
};
#endif
