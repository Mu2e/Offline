#ifndef TStnInputModule_hh
#define TStnInputModule_hh

#include "TList.h"
#include "TStnModule.hh"
#include "TList.h"

class TStnEvent;
class TStnNode;
class TChain;
class TFile;
class TStnDataset;

class TStnInputModule: public TStnModule {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TFile*                fFile;
  bool                  fOwnChain;      // Is it our chain?
  TChain*               fChain;		// ! pointer to the analyzed TChain
  Int_t                 fCurrent;	// current Tree number in a TChain
  Double_t              fEntry;		// entry number in the chain
  TList*                fDatasetList;   // list of datasets to be processed
  Int_t                 fSplitInd;	// split number of fSplitTot
  Int_t                 fSplitTot;	// total number of splits
  Int_t                 fNEntries;      // number of events on input
  Int_t                 fNFiles;        // number of filess on input
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnInputModule(const char* name  = "InputModule", 
		  const char* title = "Input Module");

  virtual ~TStnInputModule();

  virtual int       BeginJob    ();
  virtual int       BeginRun    ();
  virtual int       Event       (Int_t i);
  virtual int       EndRun      ();
  virtual int       EndJob      ();

  int               InitChain(const char* FileName, const char* TreeName);
  virtual int       AddDataset(TStnDataset* Dataset, int Print = 0);
  virtual int       RegisterInputBranches(TStnEvent* Event) = 0;
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  virtual Double_t  GetEntries () { return fNEntries; }
  virtual Double_t  GetFiles   () { return fNFiles; }
  Double_t          GetEntry   () { return fEntry; }
  Int_t             FindEvent  (Int_t Run, Int_t Event);
  TStnDataset*      GetDataset (int i) {
    return (TStnDataset*) fDatasetList->At(i); }
  TChain*           GetChain   () { return fChain; }

					// returns pointer to TStnDataBlock,
					// but don't want to do type casting

  virtual TStnNode* GetNode(const char* BranchName,  const char* ClassName);
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual Int_t     SetBranches() { return 0; }
  // this jobs is the ind'th job of tot jobs, split data accordingly
  virtual void      SetSplit(Int_t ind, Int_t tot) {
    fSplitInd = ind;
    fSplitTot = tot;
  }

//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual Int_t     NextEvent(Int_t IEntry);
  virtual Int_t     LoadEntry(Int_t IEntry);

  ClassDef(TStnInputModule,0)   // Base class for STNTUPLE input module
};
#endif
