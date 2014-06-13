#ifndef TStnDataset_hh
#define TStnDataset_hh

#include "TNamed.h"
#include "TString.h"
#include "TObjArray.h"

class   TChain;
class   TCdf2Files;

class TStnDataset: public TNamed {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TChain*     fChain;
  TObjArray*  fListOfFilesets;	// metadata: list of filesets
  TObjArray*  fListOfFiles;	// metadata: list of files
  TObjArray*  fListOfBadFiles;  // list of bad file names
  Double_t    fNEvents;         // events in the chain, bad files not counted
  Int_t       fPrintLevel;
  Int_t       fMinRunNumber;
  Int_t       fMaxRunNumber;
  TString     fBook;
  Int_t       fCataloged;	// flag: 1 is metadata available, 0 otherwise
  Int_t       fMcFlag;
  TString     fFilesetFormat;
  Int_t       fDoneBadFiles;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// non-cataloged dataset
  TStnDataset(const char* Name = "");
					// cataloged dataset
  TStnDataset(const char* Book  ,
	      const char* Name  ,
	      Int_t       MinRun = 1, 
	      Int_t       MaxRun = 99999999,      // hopefully we wont reach
	      const char* Type   = "STNTUPLE");

  TStnDataset(const char* Book                ,
	      const char* Name                ,
	      const char* Fileset             ,
	      const char* File    = ""        ,
	      const char* Type    = "STNTUPLE");

  virtual ~TStnDataset();

  Int_t Init(const char* Book  ,
	     const char* Name  ,
	     Int_t       MinRun = 1, 
	     Int_t       MaxRun = 99999999,
	     const char* Type  = "STNTUPLE");
//-----------------------------------------------------------------------------
// dataset = "book:dataset:fileset:file"
// accessors
//-----------------------------------------------------------------------------
  Int_t       GetPrintLevel   () { return fPrintLevel; }
  TChain*     GetChain        () { return fChain;      }
  Int_t       GetNFilesets    () { return fListOfFilesets->GetEntriesFast(); }
  Int_t       GetNFiles       () { return fListOfFiles->GetEntriesFast   (); }
  Double_t    GetNEvents      () { return fNEvents;    }
  Int_t       GetCataloged    () { return fCataloged;  }
  Int_t       GetMinRunNumber () { return fMinRunNumber; }
  Int_t       GetMaxRunNumber () { return fMaxRunNumber; }
  const char* GetBook         () { return fBook.Data() ; }
  Int_t       GetMcFlag       () { return fMcFlag;       }
  const char* GetFilesetFormat() { return fFilesetFormat.Data(); }

				// this is list of TCdf2Files structures

  TObjArray*  GetListOfFilesets() { return fListOfFilesets; }
  TObjArray*  GetListOfFiles   () { return fListOfFiles;    }
  TObjArray*  GetListOfBadFiles() { return fListOfBadFiles; }
  Int_t       DoneBadFiles     () { return fDoneBadFiles;   }

  TCdf2Files* GetGdf2Files(Int_t I) { 
    return (TCdf2Files*) fListOfFiles->UncheckedAt(I); 
  }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void  SetPrintLevel(Int_t Level) { fPrintLevel = Level; }
  void  SetCataloged (Int_t Flag ) { fCataloged  = Flag ; }
  void  SetMcFlag    (Int_t Flag ) { fMcFlag     = Flag ; }
  void  SetDoneBadFiles  (Int_t Flag = 1) { fDoneBadFiles = Flag;   }
//-----------------------------------------------------------------------------
// add file to non-cataloged dataset, file name - fully specified
//-----------------------------------------------------------------------------
  Int_t AddFile      (const char* Name);
//-----------------------------------------------------------------------------
// add file to a cataloged dataset, assume list of bad files to be read 
// in from the very beginning
//-----------------------------------------------------------------------------
  Int_t AddFile      (const char* Name,  const char* FilesetName,
		      Float_t     Size,  Int_t       NEvents    ,
		      Int_t       LoEvt, Int_t       LoRun      ,
		      Int_t       HiEvt, Int_t       HiRun      ,
		      Int_t       StatusCode = 0);

  Int_t AddFileset   (const char* Name);

  // skip this file when chaining the dataset
  Int_t AddBadFile   (const char* Name);

//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Clear(Option_t* Opt = "");
  void  Print(Option_t* Opt = "") const ;

  ClassDef(TStnDataset,0)
};

#endif
