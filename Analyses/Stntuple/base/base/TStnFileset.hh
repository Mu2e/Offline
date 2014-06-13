#ifndef TStnFileset_hh
#define TStnFileset_hh

#include "TNamed.h"
#include "TString.h"
#include "TObjArray.h"
#include "TUrl.h"

class   TCdf2Files;

class TStnFileset: public TNamed {
public:
  enum {
    kDisk     = 0,		// fileset stored on non-DCache data server
    kDCache   = 1		// fileset stored in DCache
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TUrl*       fDataServer;	// name of the host where the fileset is stored
  Int_t       fLocation;
  Int_t       fNFiles;
  Int_t       fNEvents;
  Int_t       fMinRunNumber;
  Int_t       fMaxRunNumber;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TStnFileset(const char* Name     = "UNINITITLISED", 
	      const char* Server   = "root://nowhere.fnal.gov//empty/dir",
	      Int_t       Location = 0,
	      Int_t       NEvents  = -1,
	      Int_t       MinRun   = -1     , 
	      Int_t       MaxRun   = 10000000);

  virtual ~TStnFileset();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TUrl*     GetDataServer  () { return  fDataServer;   }
  Int_t     GetLocation    () { return  fLocation;     }
  Int_t     GetNFiles      () { return  fNFiles;       }
  Double_t  GetNEvents     () { return  fNEvents;      }
  Double_t  GetMinRunNumber() { return  fMinRunNumber; }
  Double_t  GetMaxRunNumber() { return  fMaxRunNumber; }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void Set(const char* Server, Int_t Location, Int_t NEvents, Int_t MinRun, 
	   Int_t MaxRun);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Clear(Option_t* Opt = "");
  void  Print(Option_t* Opt = "") const ;

  ClassDef(TStnFileset,0)
};

#endif
