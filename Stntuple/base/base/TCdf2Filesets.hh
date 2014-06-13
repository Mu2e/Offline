///////////////////////////////////////////////////////////////////////////////
// root [14] db.DescribeTable("cdf2_filesets")
// FILESET_NAME                   VARCHAR2
// CREATE_TIME                    NUMBER
// DS_NAME_ID                     VARCHAR2
// TAPE_LABEL                     VARCHAR2
// FILE_COUNT                     NUMBER
// TAPE_PARTITION                 NUMBER
///////////////////////////////////////////////////////////////////////////////
#ifndef TCDF2_FILESETS_HH
#define TCDF2_FILESETS_HH

#include "TObject.h"
#include "TString.h"

class TCdf2Filesets: public TObject {
public:
  TString  fFILESET_NAME;
  Int_t    fCREATE_TIME;
  TString  fDS_NAME_ID;
  TString  fTAPE_LABEL;
  Int_t    fFILE_COUNT;
  Int_t    fTAPE_PARTITION;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TCdf2Filesets();
  ~TCdf2Filesets();
				// ****** overloaded methods of TObject

  void Clear(Option_t* Opt="");
  void Print(Option_t* Opt="") const;

  ClassDef(TCdf2Filesets,1)
};
#endif
