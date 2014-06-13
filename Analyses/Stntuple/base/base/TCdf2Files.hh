///////////////////////////////////////////////////////////////////////////////
// root [14] db.DescribeTable("cdf2_files")
// DS_NAME_ID                     VARCHAR2
// FILE_NAME                      VARCHAR2
// CREATION_TIME                  NUMBER
// FILE_SIZE                      NUMBER
// EVENT_COUNT                    NUMBER
// LOW_EVENT                      NUMBER
// LOW_RUN                        NUMBER
// HIGH_EVENT                     NUMBER
// HIGH_RUN                       NUMBER
// LUM_SUM_ONLINE                 NUMBER
// LUM_SUM_OFFLINE                NUMBER
// CREATE_USER                    VARCHAR2
// CREATE_DATE                    NUMBER
// UPDATE_USER                    VARCHAR2
// UPDATE_DATE                    NUMBER
// FILESET_NAME                   VARCHAR2
// STATUS                         NUMBER
///////////////////////////////////////////////////////////////////////////////
#ifndef TCDF2_FILES_HH
#define TCDF2_FILES_HH

#include "TObject.h"
#include "TString.h"

class TCdf2Files: public TObject {
public:
  TString   fDS_NAME_ID;
  TString   fFILE_NAME;
  Int_t     fCREATION_TIME;
  Float_t   fFILE_SIZE;
  Int_t     fEVENT_COUNT;
  Int_t     fLOW_EVENT;
  Int_t     fLOW_RUN;
  Int_t     fHIGH_EVENT;
  Int_t     fHIGH_RUN;
  Float_t   fLUM_SUM_ONLINE;
  Float_t   fLUM_SUM_OFFLINE;
  TString   fCREATE_USER;
  Int_t     fCREATE_DATE;
  TString   fUPDATE_USER;
  Int_t     fUPDATE_DATE;
  TString   fFILESET_NAME;
  Int_t     fSTATUS;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TCdf2Files();
  ~TCdf2Files();
				// ****** overloaded methods of TObject

  void Clear(Option_t* Opt="");
  void Print(Option_t* Opt="") const;

  ClassDef(TCdf2Files,1)
};
#endif
