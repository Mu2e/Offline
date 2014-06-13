#ifndef STNTUPLE_TStnModule_hh
#define STNTUPLE_TStnModule_hh

#include "TNamed.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFolder.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

class TStnAna;
class TCanvas;
class TStnHeaderBlock;
class TStnDataBlock;
class TStnNode;
class TStnEvent;
class TStnGoodRunList;

class TStnModule: public TNamed {
public:

  enum {
    kDisabled = 0,
    kFilter   = 1,
    kVeto     = 2
  };

  enum { kNDebugBits = 100 };

protected:
  int              fEnabled;
  int              fInitialized;
  int              fLastRun;
  Int_t            fFilteringMode;	   // 1 if module is used as a filter
  Int_t            fPassed;		   // 1 if event passed processing
  Int_t            fPrintLevel;		   // print or debug level
  Int_t            fMyronFlag;
  TStnAna*         fAna;		   // ! backward pointer to TStnAna
  TFolder*         fFolder;		   // ! owned by the module, don't write
  TObjArray*       fListOfL3TrigNames;     // ! list of L3 trigger names
  TObjArray*       fListOfL3Triggers;      // ! list of passed L3 triggers
  int              fDebugBit[kNDebugBits]; // ! hopefully, it will be enough
public:
  TStnModule();
  TStnModule(const char* name, const char* title);
  virtual ~TStnModule();

  virtual int BeginJob    ();
  virtual int BeginRun    ();
  virtual int Event       (Int_t i);
  virtual int EndRun      ();
  virtual int EndJob      ();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int              GetInitialized     () { return fInitialized;   }
  Int_t            GetFilteringMode   () { return fFilteringMode; }
  int              GetEnabled         () { return fEnabled;       }
  int              GetLastRun         () { return fLastRun;       }
  int              GetDebugBit   (int I) { return fDebugBit[I];   }

  TObjArray*       GetListOfHistograms() { 
    Warning("GetListOfHistograms",Form(" from %s\n",GetName()));
    Warning("GetListOfHistograms",
 	    "2003.02.13: obsolete, use TStnModule::DeleteHistograms(). Thanks, Pasha");
    DeleteHistograms(fFolder);
    return new TObjArray;
  }

  TObjArray*       GetListOfL3TrigNames() { return fListOfL3TrigNames; }
  TObjArray*       GetListOfL3Triggers () { return fListOfL3Triggers;  }

  Int_t            PrintLevel         () { return fPrintLevel;    }
  TStnAna*         GetAna             () { return fAna;           }
  Int_t            GetMyronFlag       () { return fMyronFlag;     }
  TFolder*         GetFolder          () { return fFolder;        }

  Int_t            GetPassed       () { 
    return (((fFilteringMode == kDisabled)            ) ||
	    ((fFilteringMode == kFilter  ) &&  fPassed) || 
	    ((fFilteringMode == kVeto    ) && !fPassed)    ); 
  }
//-----------------------------------------------------------------------------
//  these methods are currently provided by TStnAna, hide them
//-----------------------------------------------------------------------------
  TStnHeaderBlock* GetHeaderBlock  ();
  TStnEvent*       GetEvent        ();
  TStnGoodRunList* GetGoodRunList  ();
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void     SetAna          (TStnAna* ana    ) { fAna           = ana;     }
  void     SetEnabled      (int      enabled) { fEnabled       = enabled; }
  void     SetFilteringMode(int      mode   ) { fFilteringMode = mode;    }
  void     SetInitialized  (int      init   ) { fInitialized   = init;    }
  void     SetLastRun      (int      run    ) { fLastRun       = run;     }
  void     SetPrintLevel   (int      level  ) { fPrintLevel    = level;   }
  void     SetMyronFlag    (int      flag   ) { fMyronFlag     = flag;    }
  void     SetPassed       (Int_t    Passed ) { fPassed        = Passed;  }
  void     SetDebugBit     (int I, int Value) { fDebugBit[I]   = Value;   }

  void    AddL3TriggerName    (const char* L3Path) {
    fListOfL3TrigNames->Add(new TObjString(L3Path)); 
  }
//-----------------------------------------------------------------------------
// the following helper methods allow to save 1 line per request, which in 
// case of 100's histograms booked is a non-negligible number
//-----------------------------------------------------------------------------
  void  DeleteHistograms(TFolder* Folder = (TFolder*) -1);

  void  AddHistogram(TObject* hist, const char* FolderName = "Hist");

  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook1D(TH1D*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");

  void  HBook2F(TH2F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		Int_t Ny, Double_t YMin, Double_t YMax,
		const char* FolderName = "Hist");

  void  HProf (TProfile*& Hist, const char* Name, const char* Title,
	       Int_t Nx, Double_t XMin, Double_t XMax,
	       Double_t YMin, Double_t YMax,
	       const char* FolderName = "Hist");

					// returns TStnDataBlock*

  void*  RegisterDataBlock(const char* BranchName, const char* ClassName);
  void   RegisterDataBlock(const char* BranchName, void* DataBlock);

  void   RegisterDataBlock(const char* BranchName, 
			   const char* ClassName, 
			   void*       DataBlock);

  Int_t  AddOutputBlock   (const char* BranchName, TStnDataBlock* Block);


  TCanvas* NewSlide(const char* name, const char* title, int nx, int ny);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Clear (const char* Opt);
  virtual void  Delete(const char* Opt);
  virtual void  Print (const char* Opt) const ;

  ClassDef(TStnModule,0)   // Base Class for the STNTUPLE analysis module
};

#endif
