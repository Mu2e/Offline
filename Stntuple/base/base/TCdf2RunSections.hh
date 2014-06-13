///////////////////////////////////////////////////////////////////////////////
// root [14] db.DescribeTable("cdf2_runsections")
// ID                             NUMBER
// RUN_NUMBER                     NUMBER
// SECTION_NUMBER                 NUMBER
// BIRTH_TIME                     NUMBER
// LOW_EVENT                      NUMBER
// HIGH_EVENT                     NUMBER
// LUM_AVERAGE_ONLINE             NUMBER
// LUM_INTEGRAL_ONLINE            NUMBER
// LUM_AVERAGE_OFFLINE            NUMBER
// LUM_INTEGRAL_OFFLINE           NUMBER
// DATA_QUALITY                   RAW   
// LUMINOSITY_VERSION             NUMBER
// LIVETIME                       NUMBER
// RUNTIME                        NUMBER
///////////////////////////////////////////////////////////////////////////////
#ifndef TCDF2_RUNSECTIONS_HH
#define TCDF2_RUNSECTIONS_HH

#include "TObject.h"
#include "TString.h"

class TCdf2RunSections: public TObject {
public:
  Int_t     fID;
  Int_t     fRUN_NUMBER;
  Int_t     fSECTION_NUMBER;
  Int_t     fBIRTH_TIME;
  Int_t     fLOW_EVENT;
  Int_t     fHIGH_EVENT;
  Float_t   fLUM_AVERAGE_ONLINE;
  Float_t   fLUM_INTEGRAL_ONLINE;
  Float_t   fLUM_AVERAGE_OFFLINE;
  Float_t   fLUM_INTEGRAL_OFFLINE;
  Int_t     fDATA_QUALITY;
  Int_t     fLUMINOSITY_VERSION;
  Int_t     fLIVETIME;
  Int_t     fRUNTIME;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  TCdf2RunSections();
  ~TCdf2RunSections();

  ClassDef(TCdf2RunSections,1)
};
#endif
