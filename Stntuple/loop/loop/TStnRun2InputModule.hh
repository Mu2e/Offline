#ifndef STNTUPLE_TStnRun2InputModule_hh
#define STNTUPLE_TStnRun2InputModule_hh

#include "TNamed.h"
#include "TObjArray.h"

class TStnAna;
class TCanvas;
class TStnHeaderBlock;
class TStnDataset;

#include "TStnInputModule.hh"

class TStnRun2InputModule: public TStnInputModule {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnRun2InputModule(); 
  TStnRun2InputModule(const char* FileName, const char* TreeName="STNTUPLE"); 
  TStnRun2InputModule(TChain*      Chain); 
  TStnRun2InputModule(TStnDataset* Dataset); 
//-----------------------------------------------------------------------------
// overloaded methods of TStnInputModule
//-----------------------------------------------------------------------------
  virtual ~TStnRun2InputModule();
  virtual TStnNode* GetNode(const char* BranchName, const char* ClassName);
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  Int_t       SetBranches();
//-----------------------------------------------------------------------------
// overloaded methods of TStnInputModule
//-----------------------------------------------------------------------------
  virtual Int_t     RegisterInputBranches(TStnEvent* Event);
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int BeginJob    ();
  virtual int BeginRun    ();
  virtual int Event       (Int_t i);
  virtual int EndRun      ();
  virtual int EndJob      ();

  ClassDef(TStnRun2InputModule,0)	// RUN II input module
};

#endif
