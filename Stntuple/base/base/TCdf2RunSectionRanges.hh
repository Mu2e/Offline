///////////////////////////////////////////////////////////////////////////////
// root [14] db.DescribeTable("cdf2_runsection_ranges")
// ID                             NUMBER
// FILE_NAME                      VARCHAR2
// RUN_NUMBER_LOW                 NUMBER
// SECTION_NUMBER_LOW             NUMBER
// RUN_NUMBER_HIGH                NUMBER
// SECTION_NUMBER_HIGH            NUMBER
///////////////////////////////////////////////////////////////////////////////
#ifndef TCDF2_RUNSECTION_RANGES_HH
#define TCDF2_RUNSECTION_RANGES_HH

#include "TObject.h"
#include "TString.h"

class TCdf2RunSectionRanges: public TObject {
public:
  Int_t     fID;
  TString   fFILE_NAME;
  Int_t     fRUN_NUMBER_LOW;
  Int_t     fSECTION_NUMBER_LOW;
  Int_t     fRUN_NUMBER_HIGH;
  Int_t     fSECTION_NUMBER_HIGH;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TCdf2RunSectionRanges();
  ~TCdf2RunSectionRanges();

  ClassDef(TCdf2RunSectionRanges,1)
};
#endif
