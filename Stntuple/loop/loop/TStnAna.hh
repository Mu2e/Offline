#ifndef STNTUPLE_TStnAna_hh
#define STNTUPLE_TStnAna_hh

#ifdef __GNUG__
#include <limits.h>
#else
#include <limits>
#endif

#include "TList.h"
#include "TFile.h"
#include "TChain.h"
#include "TFolder.h"
#include "TProfile.h"

class TStnNode;
class TStnEvent;
class TStnModule;
class TStnInputModule;
class TStnOutputModule;
class TStnHeaderBlock;
class TStnDataBlock;
class TStnDBManager;
class TStnGoodRunList;
class TEventList;
class TStnDataset;
class TStnRunSummary;
class TVisManager;

class TStnAna : public TNamed {

  struct EventList_t {
    Int_t  fRun;
    Int_t  fEvent;
    Int_t  fFound;
  };

//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:

  TString           fMode;              // mode = "run1", "run2" or "gen"
  Int_t             fInitialized;       //
  TList*            fModuleList;        //
  int               fRunNumber;         // 
  int               fSectionNumber;     // ! added to rev 1.37
  int               fEventNumber;       //
  Int_t             fNPassedEvents;     //
  Int_t             fNProcessedEvents;  //
  Int_t             fNEventsToReport;   //
  TStnInputModule*  fInputModule;       //
  Double_t          fEntry;		// entry number in the chain

  TFile*            fOutputFile;	// output file (NULL by default)
  TTree*            fOutputTree;	// output tree

  TFolder*          fFolder;	        // "Ana" folder
					// event always has an event header 
					// block, this block contains event/run
					// numbers
  TStnHeaderBlock*  fHeaderBlock;	// pointer to event header proxy

					// event is created by the input module
  TStnEvent*        fEvent;		// pointer to event (welcome back!)

  TStnOutputModule* fOutputModule;

  TStnDBManager*    fDBManager;

  TStnGoodRunList*  fGoodRunList;       // NULL by default
  Int_t             fGoodRun;		// current run flag
  TList*            fListOfRuns;	// list of processed GOOD runs (owned)
  Int_t             fMinRunNumber;
  Int_t             fMaxRunNumber;
  Int_t             fPrintLevel;
  Int_t             fMcFlag;            // !

  EventList_t*      fEventList;	        // ! list of events to be processed
//-----------------------------------------------------------------------------
// visualization hook
//-----------------------------------------------------------------------------
  TVisManager*      fVisManager;	// vis. manager. default - NULL

  Int_t (*fDisplayEventRoutine)(TVisManager* Vm);
  Int_t (*fSetTitleNodeRoutine)(TVisManager* Vm, TStnHeaderBlock* Block);
//-----------------------------------------------------------------------------
//  luminosity histograms
//-----------------------------------------------------------------------------
  TProfile*         fIntLumiTev;
  TProfile*         fIntLumiLive;
  TProfile*         fIntLumiOffl;
//-----------------------------------------------------------------------------
// internal memory management
//-----------------------------------------------------------------------------
  class Cleaner {
  public: 
    Cleaner();
    ~Cleaner();
  };

  friend class Cleaner;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// mode = "run1", "run2" or "gen"
  TStnAna(const char*  filename=0, const char* mode = "run2");
  TStnAna(TChain*      chain     , const char* mode = "run2");
  TStnAna(TStnDataset* Dataset   , const char* mode = "run2");

  ~TStnAna();
				// initialization
  Int_t     Init();
//-----------------------------------------------------------------------------
// by default filtering is OFF
//-----------------------------------------------------------------------------
  int         AddModule   (TStnModule* m, Int_t FilteringMode = 0);

  TStnModule* AddModule   (const char* ClassName, 
			   int         FilteringMode = 0,
			   const char* ModuleName    = 0,
			   const char* Title         = 0);

  int         DeleteModule(const char* name);
  int         ReloadModule(const char* name, const char* filename = 0);
  int         Enable      (const char* name);
  int         Disable     (const char* name);

  virtual int BeginJob        ();
  virtual int BeginRun        ();
  virtual int EndJob          ();

  virtual int Run             (Int_t Nev    =  0     ,
			       Int_t RunMin = -1     ,
			       Int_t RunMax = INT_MAX,
			       Int_t StartEntry = 0);

  virtual int Continue        (Int_t Nev   );
  virtual int ProcessRun      (Int_t RunMin, Int_t RunMax = -1);
  virtual int ProcessEntry    (Int_t Ientry);
  virtual int ProcessEvent    (Int_t Run, Int_t Event);
  virtual int ProcessEventList(TEventList*   EventList);
  virtual int ProcessEventList(Int_t*        EventList);
  virtual int SetSplit        (Int_t ind, Int_t tot); // run part ind of tot
  virtual int AddDataset      (TStnDataset* d);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TStnGoodRunList*  GetGoodRunList  () { return fGoodRunList;  }
  Int_t             GetInitialized  () { return fInitialized;  }
  TStnInputModule*  GetInputModule  () { return fInputModule;  }
  TStnOutputModule* GetOutputModule () { return fOutputModule; }
  TStnHeaderBlock*  GetHeaderBlock  () { return fHeaderBlock;  }
  TList*            GetListOfModules() { return fModuleList;   }
  TStnDBManager*    GetDBManager    () { return fDBManager;    }
  TVisManager*      GetVisManager   () { return fVisManager;   }

  Int_t             NProcessedEvents() { return fNProcessedEvents; }
  Int_t             NPassedEvents   () { return fNPassedEvents;    }

  TStnModule* GetModule   (const char* name) {
    return (TStnModule*) fModuleList->FindObject(name);
  }

  void*  RegisterDataBlock(const char* BranchName, const char* ClassName);

				// users: DataBlock should be of TStnDataBlock**

  void  RegisterDataBlock(const char* BranchName, void* DataBlock);
//-----------------------------------------------------------------------------
//  finally: this is the right signature!
//-----------------------------------------------------------------------------
  void  RegisterDataBlock(const char* BranchName, 
			  const char* ClassName, 
			  void*       DataBlock);

  Double_t   GetEntry           () { return fEntry;           }
  TStnEvent* GetEvent           () { return fEvent;           }
  TFolder*   GetFolder          () { return fFolder;          }
  TList*     GetModuleList      () { return fModuleList;      }
  TList*     GetListOfRuns      () { return fListOfRuns;      }
  Int_t      GetNEventsToReport () { return fNEventsToReport; }
  Int_t      GetMcFlag          () { return fMcFlag;          }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetInputModule    (TStnInputModule*  M  );
  void  SetOutputModule   (TStnOutputModule* M  );
  void  SetGoodRunList    (TStnGoodRunList*  Grl) { fGoodRunList     = Grl; }
  void  SetEvent          (TStnEvent*        Ev ) { fEvent           = Ev;  }
  void  SetDBManager      (TStnDBManager*    Dbm) { fDBManager       = Dbm; }
  void  SetVisManager     (TVisManager*      Vm ) { fVisManager      = Vm;  }
  void  SetNEventsToReport(Int_t             N  ) { fNEventsToReport = N;   }
  void  SetPrintLevel     (Int_t             L  ) { fPrintLevel      = L;   }
  Int_t SetOutputFile     (const char* Filename );
  void  SetEventList      (Int_t*      EventList);
//-----------------------------------------------------------------------------
// set callback routines
//-----------------------------------------------------------------------------
  int SetDisplayEventRoutine(Int_t (*f)(TVisManager*)) {
    fDisplayEventRoutine = f;
    return 0;
  }

  int SetSetTitleNodeRoutine(Int_t (*f)(TVisManager*, TStnHeaderBlock*)) {
    fSetTitleNodeRoutine = f;
    return 0;
  }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void          Help       (const char* Item  = 0);
  void          PrintStat  (Int_t NEvents, const char* BranchName = "");
  int           SaveFolder (TFolder* Folder, TDirectory* Dir);
  void          SaveHist   (const char* Filename, Int_t Mode = 2);
  Int_t         MergeHistograms(const char* ListOfFiles, 
				const char* OutputFile);
//-----------------------------------------------------------------------------
// visualization
//-----------------------------------------------------------------------------
  int DisplayEvent() {
    if (fDisplayEventRoutine && fVisManager) { 
      return fDisplayEventRoutine(fVisManager);
    }
    else return 0;
  }

  int SetTitleNode() {
    if (fSetTitleNodeRoutine && fVisManager) {
      return fSetTitleNodeRoutine(fVisManager,fHeaderBlock);
    }
    else return 0;
  }
//-----------------------------------------------------------------------------
//  overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Clear (const char* Opt = "");
  void  Print (const char* Opt = "") const ;
  
protected:

  Int_t  NBytesRead(TBranch* Branch, Double_t& TotBytes, Double_t& ZipBytes);
  Int_t  AddFolders(TFolder*   Fol1, TFolder*   Fol2);
  Int_t  AddArrays (TObjArray* A1  , TObjArray* A2  );

  ClassDef(TStnAna,0)  // STNTUPLE event loop utility
};
#endif
