//_____________________________________________________________________________
//BEGIN_HTML
// TStnAna provides a simple utility to loop over the events stored in 
// STNTUPLE. See *_ana.C scripts in <b><a href=
// "http://purdue-cdf.fnal.gov/CdfCode/source/Stntuple/ana/">Stntuple/ana</a>
// </b> directory for the examples of its use
//END_HTML
//
//  fPrintLevel = 10: print lumi statistics on all the runs
//  fPrintLevel = 11: print lumi statistics on good runs only
//_____________________________________________________________________________
#include "cstdlib"
#include <cmath>

#include "TROOT.h"
#include "TChain.h"
#include "TInterpreter.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TLeafObject.h"
#include "TObjString.h"
#include "TFolder.h"
#include "TH1.h"
#include "TEventList.h"

#include "Stntuple/base/TStnDataset.hh"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnDBManager.hh"
#include "Stntuple/obj/TStnRunSummary.hh"
#include "Stntuple/obj/TStnGoodRunList.hh"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/loop/TStnRun2InputModule.hh"
#ifndef STNTUPLE_NOSAM
#include "Stntuple/loop/TStnSamInputModule.hh"
#endif
#include "Stntuple/loop/TStnOutputModule.hh"

ClassImp(TStnAna)
//_____________________________________________________________________________
TStnAna::TStnAna(const char* Filename, const char* Mode) : 
  TNamed ("StnAna","STNTUPLE event loop utility"),
  fMode  (Mode)
{
  // Input parameters:
  // -----------------
  // filename: a name of the ROOT (.root) file with the STNTUPLE tree
  //           if `filename' is a directory, then all the .root files are 
  //           chained
  // Mode    : "run2" - working with RUN2 STNTUPLE's
  //           the rest - not implemented yet

  Init();
  

  fInputModule = NULL;
  if      (fMode == "run2") {
    TString fn(Filename);
    if( fn.CompareTo("sam",TString::kIgnoreCase)==0  ) {
#ifdef STNTUPLE_NOSAM
      printf("ERROR: SAM requested, but SAM input module is not linked.\n");
#else
      fInputModule = new TStnSamInputModule();
#endif
    } else {
      fInputModule = new TStnRun2InputModule(Filename);
    }
    fInputModule->SetAna(this);
  }

}

//_____________________________________________________________________________
TStnAna::TStnAna(TChain* Chain, const char* Mode)  : 
  TNamed ("StnAna","STNTUPLE event loop utility"),
  fMode(Mode)
{
  // Input parameters: 
  // -----------------
  // Chain : input chain, as it is created outside, it also has to be deleted
  //         outside

  if (! Chain) {
    Error("TStnAna","input tree is not defined");
    return;
  }

  Init();

  if      (fMode == "run2") {
    fInputModule = new TStnRun2InputModule(Chain);
    fInputModule->SetAna(this);
  }
}

//_____________________________________________________________________________
TStnAna::TStnAna(TStnDataset* Dataset, const char* Mode)  : 
  TNamed ("StnAna","STNTUPLE event loop utility"),
  fMode  (Mode)
{
  // Input parameters: 
  // -----------------
  // Dataset : input dataset created outside. As such it also has to be deleted
  //           outside

  if (! Dataset) {
    Error("TStnAna","input tree is not defined");
    return;
  }

  Init();

  if      (fMode == "run2") {
    fInputModule = new TStnRun2InputModule(Dataset);
    fInputModule->SetAna(this);
  }

  fMcFlag = Dataset->GetMcFlag();
}

//_____________________________________________________________________________
int TStnAna::Init() {
  fInitialized         = 0;
  fModuleList          = new TList();
  fEventNumber         = -1;
  fSectionNumber       = -1;
  fRunNumber           = -1;
  fMinRunNumber        = -1;
  fMaxRunNumber        = INT_MAX;
  fEvent               = 0;
  fNProcessedEvents    = 0;
  fNPassedEvents       = 0;
  fNEventsToReport     = 500;
  fOutputFile          = 0;
  fOutputTree          = 0;
  fOutputModule        = 0;
  fDBManager           = TStnDBManager::Instance();
  fVisManager          = 0;
  fDisplayEventRoutine = 0;
  fSetTitleNodeRoutine = 0;
  
  TObject* o = gROOT->GetRootFolder()->FindObject("Ana");
  if (o != 0) {
    gROOT->GetRootFolder()->Remove(o);
    o->Delete();
  }

  fFolder = gROOT->GetRootFolder()->AddFolder("Ana","TStnAna Folder");

  TH1::AddDirectory(0);

  //  fListOfGoodRuns = 0;

  fGoodRun        = 1;
  fListOfRuns     = new TList();
  fPrintLevel     = 0;
  fGoodRunList    = 0;
  fMcFlag         = 0;

  fEventList      = 0;

  return 0;
}


//_____________________________________________________________________________
TStnAna::~TStnAna()
{
  // Get ride of the event
  //printf(" TStnAna: event.\n"); fflush(stdout); fflush(stderr);
  delete fEvent;

  fListOfRuns->Delete();
  delete fListOfRuns;

//-----------------------------------------------------------------------------
// so far deleting a folder leads to deletion of the histograms as well...
//-----------------------------------------------------------------------------
  //printf(" TStnAna: folder.\n"); fflush(stdout); fflush(stderr);
  gROOT->GetRootFolder()->Remove(fFolder);
  fFolder->Clear();
  delete fFolder;

  //printf(" TStnAna: mods.\n"); fflush(stdout); fflush(stderr);
  fModuleList->Delete();
  delete fModuleList;

  //printf(" TStnAna: input.\n"); fflush(stdout); fflush(stderr);
  delete fInputModule;
  //printf(" TStnAna: output.\n"); fflush(stdout); fflush(stderr);
  delete fOutputModule;
  //printf(" TStnAna: output file.\n"); fflush(stdout); fflush(stderr);
  delete fOutputFile;

  //printf(" TStnAna: eventlist.\n"); fflush(stdout); fflush(stderr);
  delete fEventList;


  // TStnAna doesn't create the good run list, it is not its job to delete it
  //printf(" TStnAna: goodruns.\n"); fflush(stdout); fflush(stderr);
//   if (fGoodRunList)
//     delete fGoodRunList;

  // Get rid of the profile histograms
  //printf(" TStnAna: profiles.\n"); fflush(stdout); fflush(stderr);
  delete fIntLumiTev;
  delete fIntLumiLive;
  delete fIntLumiOffl;
}


//_____________________________________________________________________________
Int_t TStnAna::SetSplit(Int_t ind, Int_t tot)
{
  if(TStnInputModule* inp = GetInputModule()){
    inp->SetSplit(ind,tot);
  } else {
    Error("TStnAna","Could not find input module");
  }
}

//_____________________________________________________________________________
Int_t TStnAna::AddDataset(TStnDataset* d)
{
  if(TStnInputModule* inp = GetInputModule()){
    inp->AddDataset(d);
  } else {
    Error("TStnAna","Could not find input module");
  }
}

//_____________________________________________________________________________
Int_t TStnAna::SetOutputFile(const char* filename)
{
  if (! filename)
    return -1;

  fOutputFile = new TFile(filename,"RECREATE");

  if (! fOutputFile->IsOpen())
    return -2;

  return 0;
}


//_____________________________________________________________________________
void TStnAna::SetInputModule(TStnInputModule* module) 
{
					// module is specified explicitly, 
					// non-default name
  module->SetAna(this);
  fInputModule = module; 
}

//_____________________________________________________________________________
void TStnAna::SetOutputModule(TStnOutputModule* module) 
{
					// module is specified explicitly, 
					// non-default name
  module->SetAna(this);
  fOutputModule = module; 
}

//_____________________________________________________________________________
int TStnAna::AddModule(TStnModule* Module, Int_t FilteringMode) {

					// module is specified explicitly, 
					// non-default name
  Module->SetAna(this);
  Module->SetFilteringMode(FilteringMode);
  fModuleList->Add(Module);
  return 0;
}

//_____________________________________________________________________________
TStnModule* TStnAna::AddModule(const char* ClName    ,
			       int         FilteringMode,
			       const char* ModuleName,
			       const char* Title     )
{
  // assume that "Name" is the name of the module class and that this class
  // is already known to the interpreter. Module is created using default 
  // constructor, so 
  // newly created module is owned by *this

  TStnModule* m = 0;

  TClass* cl = gROOT->GetClass(ClName);
  if (! cl) {
    Error("AddModule","class %s is not known to the interpreter",ClName);
  }
  else {

    m = (TStnModule*) cl->New();
    if (ModuleName) m->SetName(ModuleName);
    if (Title     ) m->SetTitle(Title);

    m->SetAna(this);
    m->SetFilteringMode(FilteringMode);
    fModuleList->Add(m);
  }
  return m;
}

//_____________________________________________________________________________
int TStnAna::DeleteModule(const char* name) {
  TObject* m = fModuleList->FindObject(name);
  if (m) {
    fModuleList->Remove(m);
    delete m;
  }
  return 0;
}

//_____________________________________________________________________________
int TStnAna::ReloadModule(const char* name, const char* filename) {
  char fname[200], soname[200], cmd[200];
  char class_name[100];

  if (this) {
    printf(">>> TStnAna::ReloadModule not implemented yet\n");
    return -1;
  }

  DeleteModule(name);
					// define .so and .cc files

  sprintf(class_name,"T%sModule",name);
  sprintf(soname,"%s/./%s.so",getenv("PWD"),class_name);
  sprintf(fname ,"%s/./%s.cc",getenv("PWD"),class_name);

					// unload the library
  sprintf(cmd,".U %s",soname);
  gInterpreter->ProcessLine(cmd);
					// load modified library
  sprintf(cmd,".L %s+",fname);
  gInterpreter->ProcessLine(cmd);
					// and finally add new module

  TClass* cl = gROOT->GetClass(class_name);
//    if (! cl) {
//      printf(">>> TStnAna::ReloadModule: class %s doesn exist\n",class_name);
//      return -2;
//    }
//    else if (! cl->InheritsFrom("TStnModule")) {
//      printf(">>> TStnAna::ReloadModule: class %s doesn inherit from TStnModule",
//  	   class_name);
//      return -3;
//    }

  AddModule((TStnModule*) cl->New());

  return 0;
}

//_____________________________________________________________________________
// is called when teh HeaderBlock is read in and the run number is not equal 
// to the last run number read
//-----------------------------------------------------------------------------
int TStnAna::BeginRun() {

  fRunNumber     = fHeaderBlock->RunNumber();
  fSectionNumber = fHeaderBlock->SectionNumber();
//-----------------------------------------------------------------------------
// initialize DB manager
//-----------------------------------------------------------------------------
  TFile*      file;
  TDirectory* olddir;
  TChain*     chain;

  chain  = fInputModule->GetChain();
//-----------------------------------------------------------------------------
//  MC generator module returns chain=0 
//-----------------------------------------------------------------------------
  if (chain) {
    file   = chain->GetFile();
    olddir = gDirectory;

    if (file && file->cd("db")) {
				// DB directory exists, read the calib const

      gDirectory->cd(Form("run_%i",fRunNumber));
      fDBManager->Read("");
    }
    olddir->cd();
  }

  TStnRunSummary* new_rs = (TStnRunSummary*) fDBManager->GetTable("RunSummary");
//-----------------------------------------------------------------------------
// if use bad run list and the run is marked as bad, return...
//-----------------------------------------------------------------------------
  if (fGoodRunList != 0)
    fGoodRun = fGoodRunList->GoodRun(fRunNumber,fHeaderBlock->SectionNumber());
  else
    fGoodRun = 1;

  if (! fGoodRun)                                         return -1;
//-----------------------------------------------------------------------------
// cache new run to calculate luminosity in the end of job
//-----------------------------------------------------------------------------
  TIter it(fListOfRuns);
  int done = 0;
  TStnRunSummary* rs;

  while ((rs = (TStnRunSummary*) it.Next())) {
    if (rs->RunNumber() == new_rs->RunNumber()) {
      done = 1;
      break;
    }
    else if (rs->RunNumber() > new_rs->RunNumber()) {
      break;
    }
  }

  if (! done) {
    TStnRunSummary* run_sum = new TStnRunSummary(*new_rs);
    if (rs)
      fListOfRuns->AddBefore(rs,run_sum);
    else
      fListOfRuns->Add(run_sum);
  }
  return 0;
}


//_____________________________________________________________________________
int TStnAna::ProcessEventList(TEventList* EventList) {
  // process event list

  int   n, ientry,rc;

  n = EventList->GetN();

  for (int i=0; i<n; i++) {
    ientry = EventList->GetEntry(i);
    rc = ProcessEntry(ientry);
    if(rc!=0 && rc!=-2 && rc!=-3) return rc;
  }

  return 0;
}

//_____________________________________________________________________________
int TStnAna::ProcessEventList(Int_t* EventList) {
  // process event list 

  SetEventList(EventList);
  int rc = Run();

  return rc;
}

//_____________________________________________________________________________
void TStnAna::SetEventList(Int_t* EventList) {
  // set event list 

  int nev;
//-----------------------------------------------------------------------------
//  determine number of events
//-----------------------------------------------------------------------------
  nev = 0;

  for (int i=0; EventList[i] > 0; i+= 2) { nev++; }

  fEventList = new EventList_t[nev+1];

  for (int i=0; i<nev; i++) { 
    fEventList[i].fRun   = EventList[2*i  ];
    fEventList[i].fEvent = EventList[2*i+1];
    fEventList[i].fFound = 0;
  }

  fEventList[nev].fRun = -1;

}

//_____________________________________________________________________________
int TStnAna::ProcessEntry(int Entry) {
  // process one event - `Entry' (!!!) in the chain
  // normally RunMin = -1, in this case process all the events
  // if RunMin > 0 and RunMax = -1, process only RunMin
  // if RunMin > 0 and RunMax > 0,  process events from RunMin to RunMax 
  //                                (including RunMin and RunMax)

  int        passed, rn, rsn, ev, rc;
				// prepare to read next event (run2 case: load
				// the tree, run1 case: also read the event)

  if (! fInitialized) {
    rc = BeginJob();
    if(rc!=0) return rc;
  }

  fEntry = Entry;
  Int_t tree_entry = fInputModule->NextEvent(int(fEntry));

  if (tree_entry < 0) {
//-----------------------------------------------------------------------------
// end of the chain...
//-----------------------------------------------------------------------------
    Warning("ProcessEntry","End of chain reached.");
    return -1;
  }
//-----------------------------------------------------------------------------
// clear all the data from the previous event
//-----------------------------------------------------------------------------
  TIter nd(fEvent->GetListOfNodes());
  while (TStnNode* node = (TStnNode*) nd.Next()) {
    node->GetDataBlock()->Clear();
  }
				// always read event header branch

  rc = fHeaderBlock->GetEntry(tree_entry);
  if (rc == 0 || rc == -1) {    // non existing entry or io error
    return -4;
  }

  rn  = fHeaderBlock->RunNumber  ();
  rsn = fHeaderBlock->SectionNumber();
  ev  = fHeaderBlock->EventNumber();
  fHeaderBlock->GetEvent()->SetEventNumber(rn,ev,rsn);
  passed = 1;
//-----------------------------------------------------------------------------
//  if event list is specified
//-----------------------------------------------------------------------------
  if (fEventList != 0) {
    passed = 0;
    for (int i=0; fEventList[i].fRun > 0; i++) {
      if ((fEventList[i].fRun == rn) && (fEventList[i].fEvent == ev)) {
	passed = 1;
	fEventList[i].fFound++;
	break;
      }
    }
    if (passed == 0)                                        return 0;
  }
//-----------------------------------------------------------------------------
//  check run number and return if it is outside the range
//-----------------------------------------------------------------------------
  if (( rn < fMinRunNumber) || (rn > fMaxRunNumber ))       return -2;

  if (fHeaderBlock->RunNumber() != fRunNumber) {
//-----------------------------------------------------------------------------
// BeginRun returns -1 if the run is bad
//-----------------------------------------------------------------------------
    BeginRun();
  }
//-----------------------------------------------------------------------------
// good/bad run check is performed at the run-section level
// if BeginRun was called, the check has already been performed and 
// fSectionNumber has been set there 
// use it not to do this check 2 times
//-----------------------------------------------------------------------------
  if (fSectionNumber != rsn) {
    if (fGoodRunList != 0) fGoodRun = fGoodRunList->GoodRun(rn,rsn);
    else                   fGoodRun = 1;
  }
  if (fGoodRun <= 0)                                        return -3;
//-----------------------------------------------------------------------------
// each module is supposed to talk to an input chain and request the data 
// (branches) it needs
//-----------------------------------------------------------------------------
  TListIter it(fModuleList);
  while (TStnModule* m = (TStnModule*) it.Next()) {
    if (m->GetEnabled()) {
					// make sure BeginJob was called at 
					// least once
      if (! m->GetInitialized()) {
	m->BeginJob();
	m->SetInitialized(1);
      }
					// check if the run # has changed

      if (fHeaderBlock->RunNumber() != m->GetLastRun()) {
	m->EndRun();
	m->BeginRun();
	m->SetLastRun(fHeaderBlock->RunNumber());
      }
					// now it is safe to call the 
					// event entry point
      if (fPrintLevel > 0) {
	printf(" ---------- Event entry point called for %s\n",m->GetName());
      }
//-----------------------------------------------------------------------------
// handle fatal errors detected by the modules
//-----------------------------------------------------------------------------
      rc = m->Event(tree_entry);
      if (rc < 0) {
	fHeaderBlock->Print(Form("%s rc=%i detected by %s, skip event",
				 "TStnaAna::ProcessEntry: FATAL ERROR ",
				 rc,m->GetName()));
	passed = 0;
      }
      else {
	passed = m->GetPassed();
      }
      if (! passed) break;
    }
  }
					// see if output needs to be written
  if (passed) {
    fNPassedEvents++;
    if (fOutputModule) {
      if (fOutputModule->GetEnabled()) {
	if (! fOutputModule->GetInitialized()) {
	  fOutputModule->BeginJob();
	  fOutputModule->SetInitialized(1);
	}
	if (fHeaderBlock->RunNumber() != fOutputModule->GetLastRun()) {
	  fOutputModule->EndRun();
	  fOutputModule->BeginRun();
	}

	fOutputModule->Event(tree_entry);
      }
    }
  }

  fNProcessedEvents++;

  if (fPrintLevel == 0) {
    if (fNProcessedEvents % fNEventsToReport == 0) {
      fHeaderBlock->Print(Form(" >>> TStnAna: Processed %10i events",
			       fNProcessedEvents));
    }
  }
				// visualization hook
  DisplayEvent();

  return 0;
}

//_____________________________________________________________________________
int TStnAna::ProcessEvent(int Run, int Event) {
  // process one event with the given run/event numbers. This method is not as 
  // fast as ProcessEntry, because it starts searching from the beginning of 
  // the tree
  // in this mode ignore bad run list

  TIter it(fModuleList);
				// prepare to read next event
  int found  = 0;
  Int_t  tree_entry, nb, rc;

  if (! fInitialized) {
    rc = BeginJob();
    if (rc!=0) return rc;
  }
//-----------------------------------------------------------------------------
// clear previous event data not to think about it in the modules
//-----------------------------------------------------------------------------
  TIter nd(fEvent->GetListOfNodes());

  while (TStnNode* node = (TStnNode*) nd.Next()) {
    node->GetDataBlock()->Clear();
  }
//-----------------------------------------------------------------------------
// try to find event in the chain
//-----------------------------------------------------------------------------
  fEntry     = 0;
  tree_entry = fInputModule->GetChain()->LoadTree(int(fEntry));

  do {
    tree_entry = fInputModule->NextEvent(int(fEntry));
    if (tree_entry < 0) break;
    nb         = fHeaderBlock->GetEntry(tree_entry);
    if ((fHeaderBlock->EventNumber() == Event) && 
	(fHeaderBlock->RunNumber  () == Run  )   ) {
      found = 1;
      break;
    }
    fEntry++;
  } while (nb);

  if (! found) {
    printf(" *** run %i event %i not found\n",Run,Event);
    return -1;
  }
  else {
    fHeaderBlock->Print();
  }

  fHeaderBlock->GetEvent()->SetEventNumber(fHeaderBlock->RunNumber(),
					   fHeaderBlock->EventNumber(),
					   fHeaderBlock->SectionNumber());
//-----------------------------------------------------------------------------
// each module is supposed to request data (branches) it needs
//-----------------------------------------------------------------------------
  fEntry     = fInputModule->GetEntry();

  if (fHeaderBlock->RunNumber() != fRunNumber) {
    BeginRun();
  }

  it.Reset();
  while (TStnModule* m = (TStnModule*) it.Next()) {
    if (m->GetEnabled()) {
					// make sure BeginJob was called 
					// at least once
      if (! m->GetInitialized()) {
	m->BeginJob();
	m->SetInitialized(1);
      }
					// check if run # has changed

      if (fHeaderBlock->RunNumber() != m->GetLastRun()) {
	m->EndRun();
	m->BeginRun();
	m->SetLastRun(fHeaderBlock->RunNumber());
      }
					// and only now it is safe to call 
					// event entry point
      m->Event(tree_entry);
    }
  }
  
  if (fOutputModule) {
    if (fOutputModule->GetEnabled()) {
      if (! fOutputModule->GetInitialized()) {
	fOutputModule->BeginJob();
	fOutputModule->SetInitialized(1);
      }
      if (fHeaderBlock->RunNumber() != fOutputModule->GetLastRun()) {
	fOutputModule->EndRun();
	fOutputModule->BeginRun();
      }
      
      fOutputModule->Event(tree_entry);
    }
  }

  fNPassedEvents++;
  fNProcessedEvents++;

  if (fPrintLevel == 0) {
    if ((fNEventsToReport > 0) && (fNProcessedEvents % fNEventsToReport == 0)) {
      fHeaderBlock->Print(Form(" >>> TStnAna: Processed %10i events",
			       fNProcessedEvents));
    }
  }
				// visualization hook
  DisplayEvent();

  return 0;
}



//_____________________________________________________________________________
int TStnAna::BeginJob() 
{
  // initialization actions to be done after construction

  char  name[200], title[200];
  int rc;

  rc = fInputModule->BeginJob();
  if(rc!=0) return rc;
  fInputModule->RegisterInputBranches(fEvent);
  fInputModule->SetInitialized(1);
  fFolder->Add(fInputModule->GetFolder());
//-----------------------------------------------------------------------------
// always read in the header block
//-----------------------------------------------------------------------------
  RegisterDataBlock("HeaderBlock",&fHeaderBlock);
  TIter it(fModuleList);
				// initialization
  
  // the input module may have opened a file, make sure we cd out 
  // of it before the modules declare histograms
  gROOT->cd();

  while (TStnModule* m = (TStnModule*) it.Next()) {
    if (m->GetEnabled()) {
      m->BeginJob();
      m->SetInitialized(1);
				// also pick up module's folder...
      fFolder->Add(m->GetFolder());
    }
  }

  if (fOutputModule && fOutputModule->GetEnabled()) {

				// make sure that all the input branches made it
				// into the output list and being read in

    TIter itt(fEvent->GetListOfInputNodes());
    TStnNode*   node;
    TObjString* name;
    const char* branch_name;
    TObjArray*  drop_list = fOutputModule->GetDropList();
    TObjArray*  keep_list = fOutputModule->GetKeepList();
    int ndropped = drop_list->GetEntriesFast();
    int nkept    = keep_list->GetEntriesFast();

    while (node = (TStnNode*) itt.Next()) {
      branch_name = node->GetName();

      // If the keep list has an element we drop everything and only
      // check the keep list
      if (nkept < 1) {
	for (int i=0; i<ndropped; i++) {
	  name = (TObjString*) drop_list->UncheckedAt(i);
	  if (strcmp(name->GetString().Data(),branch_name) == 0)
	    goto NEXT_BRANCH;
	}
      }
      else {
	bool dropIt = true;
	for (int i=0; i<nkept; i++) {
	  name = (TObjString*) keep_list->UncheckedAt(i);
	  if (strcmp(name->GetString().Data(),branch_name) == 0) {
	    dropIt = false;
	    break;
	  }
	}
	if (dropIt)
	  goto NEXT_BRANCH;
      }
      
      // Add this node to the file
      if (! fEvent->GetListOfNodes()->FindObject(node))
	fEvent->GetListOfNodes()->Add(node);
      fEvent->GetListOfOutputNodes()->Add(node);

    NEXT_BRANCH:;
    }

    fOutputModule->BeginJob();
    fOutputModule->SetInitialized(1);
    fFolder->Add(fOutputModule->GetFolder());
  }
//-----------------------------------------------------------------------------
//  book luminosity histograms
//-----------------------------------------------------------------------------
  sprintf(name, "int_lumi_tev");
  sprintf(title,"Integrated delivered luminosity");
  fIntLumiTev = new TProfile(name,title,1000,110000,210000,0,1e10);

  sprintf(name, "int_lumi_tape");
  sprintf(title,"Integrated recorded luminosity");
  fIntLumiLive = new TProfile(name,title,1000,110000,210000,0,1e10);

  sprintf(name, "int_lumi_offl");
  sprintf(title,"Integrated offline luminosity");
  fIntLumiOffl = new TProfile(name,title,1000,110000,210000,0,1e10);

  fFolder->Add(fIntLumiTev);
  fFolder->Add(fIntLumiLive);
  fFolder->Add(fIntLumiOffl);
				// visualization hook
  SetTitleNode();

  fInitialized = 1;

  return 0;
}

//_____________________________________________________________________________
int TStnAna::EndJob() {
  // actions at end of job

  static int nlines = 0;
//-----------------------------------------------------------------------------
// fill luminosity histograms 
//----------------------------------------------------------------------------- 
  TIter it(fListOfRuns);

  float ltev  = 0;
  float ltape = 0;
  float loffl = 0;

  while (TStnRunSummary* rs = (TStnRunSummary*) it.Next()) {

    if (fPrintLevel == 11) {
      if (nlines == 0) {
	printf("--------------------------------------------------");
	printf("----------------------\n");
	printf("      run   rc status                   ");
	printf(" L(TeV)     L(live)   L(offline)\n");
	printf("---------------------------------------------------");
	printf("---------------------\n");
      }
      printf(" %8i  %3i %-20s  %10.3f  %10.3f   %10.3f\n",
	     rs->RunNumber(),
	     1,
	     "GOOD",
	     rs->LumiTev(),
	     rs->LumiTape(),
	     rs->OfflineLumiRS());
      nlines++;
      if (nlines == 50) nlines = 0;
    }
    ltev  += rs->LumiTev();
    ltape += rs->LumiTape();
    loffl += rs->OfflineLumiRS();
    
    fIntLumiTev->Fill (rs->RunNumber(),ltev);
    fIntLumiLive->Fill(rs->RunNumber(),ltape);
    fIntLumiOffl->Fill(rs->RunNumber(),loffl);
  }

  if (fPrintLevel >= 0) {
    printf("----- end job: ---- %s\n", GetName());
    printf(" L(TeV) : %10.3f\n",ltev);
    printf(" L(live): %10.3f\n",ltape);
    printf(" L(offl): %10.3f\n",loffl);
  }
//-----------------------------------------------------------------------------
// other actions
//----------------------------------------------------------------------------- 
  fInputModule->EndJob();

  TIter itt(fModuleList);
  while (TStnModule* m = (TStnModule*) itt.Next()) {
    if (m->GetEnabled()) {
      m->EndRun();
      m->EndJob();
    }
  }
  
  if (fOutputModule && fOutputModule->GetEnabled()) {
    fOutputModule->EndJob();
  }

  if (fPrintLevel > -2)
    printf(" >>> TStnAna::EndJob: processed %10i events, passed %10i events\n",
	   fNProcessedEvents,fNPassedEvents);

  if (fEventList) {
    printf(" >>>  strip summary:-\n");
    for (int i=0; fEventList[i].fRun>0; i++) {
      printf (" run,event : %6i,%8i   found=%i\n",
	      fEventList[i].fRun, 
	      fEventList[i].fEvent,
	      fEventList[i].fFound);
    }

    delete [] fEventList;
    fEventList = 0;
  }

  return 0;
}

//_____________________________________________________________________________
int TStnAna::Run(int NEvents, Int_t MinRunNumber, Int_t MaxRunNumber,
		 Int_t StartEntry) {
  // process NEvents events starting from the first entry in the tree/chain
  // perform termination actions (whatever they are) in the end
  int rc;
  rc = BeginJob();
  if(rc!=0) return rc;

  double nentries = fInputModule->GetEntries();

  if (NEvents > 0) {
    if (NEvents <= nentries) {
      nentries = NEvents;
    }
    else {
      Warning("Run",
	      " STNTUPLE has only %10.0f entries, those will be processed",
	      nentries);
    }
  }

  fMinRunNumber = MinRunNumber;
  fMaxRunNumber = MaxRunNumber;

  fNProcessedEvents = 0;
  fNPassedEvents    = 0;
  fEntry            = StartEntry-1;

  Continue(int(nentries));

  EndJob();

  return 0;
}


//_____________________________________________________________________________
int TStnAna::ProcessRun(int RunMin, int RunMax) {
  // process all events of a given run
  // perform termination actions (whatever they are) in the end

  int rc;

  rc = BeginJob();
  if (rc!=0) return rc;

  double nentries = fInputModule->GetEntries();

  printf(" processing events for the run range from %8i to %8i\n",
	 RunMin,RunMax);

  fNProcessedEvents = 0;
  fNPassedEvents    = 0;
  fEntry            = -1;
  fMinRunNumber     = RunMin;

  if (RunMax == -1) fMaxRunNumber = RunMin;
  else              fMaxRunNumber = RunMax;

  for (int i=0; i<nentries; i++) {
//-----------------------------------------------------------------------------
// ProcessEntry increments fEntry
//-----------------------------------------------------------------------------
    rc = ProcessEntry(i);
    if (rc!=0 && rc!=-2 && rc!=-3) break;
  }

  EndJob();

  return 0;
}


//_____________________________________________________________________________
int TStnAna::Continue(int nev) {
  // process next `nev' events w/o doing any initialization/termination
  // fEntry: last processed entry (in the chain)

  int      rc;

  Double_t i0 = fEntry+1;

  int i=0;
  int ientry = 0;

  while (1) {
//-----------------------------------------------------------------------------
// ProcessEntry increments fEntry - need a way of signalling end-of-file
// rc = -1 : problem with reading the Header, most probably - end of file
// rc = -2 : run outside the requested limits
// rc = -3 : bad run according to the good run list used
//-----------------------------------------------------------------------------
    rc = ProcessEntry(int(i0+ientry));
    if (rc!=0 && rc!=-2 && rc!=-3) {
      break;
    }
    ientry ++;
    if (rc == 0) {
      i++;
      if (i >= nev) break;
    }
  }
  return 0;
}


//_____________________________________________________________________________
void* TStnAna::RegisterDataBlock(const char*     BranchName,
				 const char*     ClassName)
{
					// make sure branch exists: the branch 
					// is a top-level branch
  TStnDataBlock* block;
  TStnNode*      node, *exists;
  TClass*        cl;

  node = (TStnNode*) fEvent->GetListOfInputNodes()->FindObject(BranchName);

  if (node) {
				// make sure we're not adding the same branch
				// twice

    exists =  (TStnNode*) fEvent->GetListOfNodes()->FindObject(BranchName);

    if (exists) node = exists;
    else fEvent->GetListOfNodes()->Add(node);

    block = node->GetDataBlock();
    block->SetNode(node);
    block->SetValid(1);
  }
  else {
    cl = gROOT->GetClass(ClassName);
    if (!cl || !cl->InheritsFrom("TStnDataBlock")) {
      Error("RegisterDataBlock",
	    "class %s doesn\'t inherit from TStnDataBlock",
	    ClassName);
      block = 0;
    }
    else {
      node  = new TStnNode("",cl,fEvent);
      fEvent->GetListOfUnusedNodes()->Add(node);
      block = node->GetDataBlock();
      block->SetNode(node);
      block->SetValid(0);
    }
  }

  return block;
}


//_____________________________________________________________________________
void TStnAna::RegisterDataBlock(const char* BranchName, void* DataBlock)
{
  // make sure branch exists: the branch is a top-level branch

  TStnDataBlock** block = (TStnDataBlock**) DataBlock;
  TStnNode        *node, *exists;

  node = (TStnNode*) fEvent->GetListOfInputNodes()->FindObject(BranchName);

  if (node) {
//-----------------------------------------------------------------------------
// make sure we're not adding the same branch twice
//-----------------------------------------------------------------------------
    exists =  (TStnNode*) fEvent->GetListOfNodes()->FindObject(BranchName);

    if (exists) node = exists;
    else fEvent->GetListOfNodes()->Add(node);

    (*block) = node->GetDataBlock();
    (*block)->SetNode(node);
    (*block)->SetValid(1);
  }
  else {
//-----------------------------------------------------------------------------
//  branch doesn't exist, no idea about the class name, 
//  so the best we can do is to create a TStnBlock...
//-----------------------------------------------------------------------------
    (*block) = new TStnDataBlock();
    (*block)->SetValid(0);
    node  = new TStnNode("",(*block),fEvent);
    (*block)->SetNode(node);
    fEvent->GetListOfUnusedNodes()->Add(node);
  }
}


//_____________________________________________________________________________
void TStnAna::RegisterDataBlock(const char* BranchName, 
				const char* ClassName, 
				void*       DataBlock)
{
  // make sure branch exists: the branch is a top-level branch

  TStnDataBlock** block = (TStnDataBlock**) DataBlock;
  TStnNode        *node, *exists;
  TClass*         cl;

  node = (TStnNode*) fEvent->GetListOfInputNodes()->FindObject(BranchName);

  if (node) {
//-----------------------------------------------------------------------------
// make sure we're not adding the same branch twice
//-----------------------------------------------------------------------------
    exists =  (TStnNode*) fEvent->GetListOfNodes()->FindObject(BranchName);

    if (exists) node = exists;
    else fEvent->GetListOfNodes()->Add(node);

    (*block) = node->GetDataBlock();
    (*block)->SetNode(node);
    (*block)->SetValid(1);
  }
  else {
//-----------------------------------------------------------------------------
//  branch doesn't exist, no idea about the class name, 
//  so the best we can do is to create a TStnBlock...
//-----------------------------------------------------------------------------
    cl = gROOT->GetClass(ClassName);
    if (!cl || !cl->InheritsFrom("TStnDataBlock")) {
      Error("RegisterDataBlock",
	    "class %s doesn\'t inherit from TStnDataBlock",
	    ClassName);
      (*block) = 0;
    }
    else {
      node  = new TStnNode("",cl,fEvent);
      fEvent->GetListOfUnusedNodes()->Add(node);
      (*block) = node->GetDataBlock();
      (*block)->SetNode(node);
      (*block)->SetValid(0);
    }
  }
}


//_____________________________________________________________________________
int TStnAna::Enable(const char* module_name) {
  TStnModule* m = (TStnModule*) fModuleList->FindObject(module_name);
  if (m) {
    m->SetEnabled(1);
    return 0;
  }
  else {
    return -1;
  }
}

//_____________________________________________________________________________
int TStnAna::Disable(const char* module_name) {
  TStnModule* m = (TStnModule*) fModuleList->FindObject(module_name);
  if (m) {
    m->SetEnabled(0);
    return 0;
  }
  else {
    return -1;
  }
}

//_____________________________________________________________________________
void TStnAna::Help(const char* item) {
  if (gSystem->Exec("which lynx >& /dev/null") == 0) {
    gSystem->Exec("xterm -geometry 100x40 lynx http://www-cdf.fnal.gov/upgrades/computing/projects/Stntuple/Stntuple.html");
  }
  else {
    printf(" *** can\'t find lynx\n");
  }
}



namespace {
				// local little struct to keep statistics
  struct BranchStat {
    //---------------------------
    // data
    //---------------------------
    TBranch*   fBranch;
    Double_t   fN;
    Double_t   fX;
    Double_t   fX2;
    Double_t   fTotBytes;
    Double_t   fZipBytes;
    Double_t   fCompRatio;
    //---------------------------
    // functions
    //---------------------------
     BranchStat(): fN(0), fX(0), fX2(0), fTotBytes(0), fZipBytes(0) {}
    //    ~BranchStat() {}

    Double_t  XMean () { return fX/fN; }
    Double_t  X2Mean() { return fX2/fN; }
    Double_t  SigXX () { return X2Mean() - XMean()*XMean(); }
    Double_t  TotBytes() {return fTotBytes;}
    Double_t  ZipBytes() {return fZipBytes;}
    Double_t  CompRatio() {return (fZipBytes>0 ? fTotBytes/fZipBytes : 0.0);}
  };
}

//_____________________________________________________________________________
Int_t TStnAna::NBytesRead(TBranch* Branch, 
			  Double_t& TotBytes, 
			  Double_t& ZipBytes) 
{
  // returns total and zipped (on disk) sizes of the `Branch'

  TotBytes += 1.0 * Branch->GetTotBytes();
  ZipBytes += 1.0 * Branch->GetZipBytes();
  TObjArray* branches = Branch->GetListOfBranches();

  int nbranches = branches->GetEntriesFast();
  if (nbranches > 0) {
    for (int i=0; i<nbranches; i++) {
      TBranch* b = (TBranch*) branches->At(i);
      NBytesRead(b,TotBytes,ZipBytes);
    }
  }
  return 0;
}

//_____________________________________________________________________________
int  TStnAna::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

  TFolder*     fol;
  TDirectory*  dir;
  TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
  Dir->cd();
  //  dir = new TDirectory(Folder->GetName(),Folder->GetName(),"");
  dir = Dir->mkdir(Folder->GetName(),Folder->GetName());
  dir->cd();

//   printf(" ------------------- Dir: %s, new dir: %s\n",
// 	 Dir->GetName(),dir->GetName());


  TIter  it(Folder->GetListOfFolders());
  while ((o = it.Next())) {
//     printf(" o->GetName, o->ClassName : %-20s %-20s\n",
// 	   o->GetName(),
// 	   o->ClassName());

    if (strcmp(o->ClassName(),"TFolder") == 0) {
      SaveFolder((TFolder*) o, dir);
      //      dir->cd();
    }
    else if (! o->InheritsFrom("TStnModule")) {
      //      printf("gDirectory->GetPath = %s\n",gDirectory->GetPath());
      o->Write();
      //      gDirectory->GetListOfKeys()->Print();
    }
  }

  Dir->cd();
  return 0;
}

//_____________________________________________________________________________
void TStnAna::SaveHist(const char* Filename, Int_t Mode) 
{
  // save histograms booked by all the modules into a file with the given name
  // Mode = 1: save folders
  // Mode = 2: save directories

  TFile* f = new TFile(Filename,"recreate");

  if (Mode == 1) {
    Error("SaveHist","Mode=1 is obsolete, use Mode=2. 2010-01-22, Pasha.");
      //    fFolder->Write();
  }
  else if (Mode == 2) {
    SaveFolder(fFolder,f);
  }

  f->Close();
  delete f;
}

//_____________________________________________________________________________
void TStnAna::PrintStat(Int_t nevents, const char* BranchName) 
{
  // Prints status information after reading `nevents' from `BranchName'
  // `BranchName = "" means all the branches
  BeginJob();
  TObjArray* list  = fInputModule->GetChain()->GetListOfBranches();

  // In case input was not properly specified
  if (! list) {
    Error("PrintStat","input chain seems empty.");
    return;
  }

  TIter it(list);
  TBranch*       branch;
  TStnNode*      node;
  TLeafObject*   leaf;

  int nbranches = list->GetEntriesFast();

  const char*       class_name;

  TStnEvent*  event = new TStnEvent();
  BranchStat* stat  = new BranchStat[nbranches+1];

  int ibb = 0;
  while (branch = (TBranch*) it.Next()) {
    if ( (strcmp(BranchName,""               ) == 0) || 
	 (strcmp(BranchName,branch->GetName()) == 0)    ) {

      if (strcmp(branch->ClassName(),"TBranchElement") == 0) {
	class_name = ((TBranchElement*)branch)->GetClassName();
      }
      else {
	leaf = (TLeafObject*) ((TBranchObject*)branch)->GetLeaf(branch->GetName());
	class_name = leaf->GetTypeName();
      }

      TClass* cl   = gROOT->GetClass(class_name);
      node         = new TStnNode(branch->GetName(),cl,event);
      branch->SetAddress(node->GetDataBlockAddress());
      branch->SetAutoDelete(kFALSE);
      event->GetListOfNodes()->Add(node);
      //    printf("%s\n",branch->GetName());
      stat[ibb].fBranch = branch;
      stat[ibb].fN      = 0;
      stat[ibb].fX      = 0;
      stat[ibb].fX2     = 0;
      ibb++;
    }
  }

  Double_t nev;
				// set the number of events to read
  nev = fInputModule->GetEntries();
  if (nev >= nevents) {
    nev = nevents;
  }
  else {
    Warning("PrintStat"," STNTUPLE has only %g entries, those will be processed",
	    nev);
  }

  Double_t nb_event, nb;
  for (int ientry=0; ientry<nev; ientry++) {
    nb_event = 0;
    event->Clear();
    for (int ib=0; ib<ibb; ib++) {
      nb = stat[ib].fBranch->GetEntry(ientry);
      stat[ib].fN  += 1;
      stat[ib].fX  += nb;
      stat[ib].fX2 += nb*nb;
      nb_event += nb;
      if (fPrintLevel == 2) {
	printf("%-48s   %10.0f bytes  \n",stat[ib].fBranch->GetName(),nb);
      }
    }
    stat[nbranches].fN  += 1;
    stat[nbranches].fX  += nb_event;
    stat[nbranches].fX2 += nb_event*nb_event;
  }

  // find compression ratio
  for (int ib=0; ib<ibb; ib++) {
    NBytesRead(stat[ib].fBranch,stat[ib].fTotBytes,stat[ib].fZipBytes);
    stat[nbranches].fTotBytes += stat[ib].TotBytes();
    stat[nbranches].fZipBytes += stat[ib].ZipBytes();
  }

  Double_t meanCompSize = 
    ( stat[nbranches].CompRatio()>0 ? 
      stat[nbranches].XMean()/stat[nbranches].CompRatio() : 0.0);

				// last piece: printing
  printf("----------------------------------------------------");
  printf("----------------------------------------------------\n");
  printf("........... branch name .....................<event size>  <sigma size> "); 
  printf(" TotBytes   ZipBytes CompFactor   %% of File\n");
  printf("----------------------------------------------------");
  printf("----------------------------------------------------\n");


  for (int ib=0; ib<ibb; ib++) {
    Double_t  diskFrac = ( stat[nbranches].ZipBytes()>0.0 ? 
	    stat[ib].ZipBytes()/stat[nbranches].ZipBytes() : 0.0 );
    printf("%-48s %7.0f %7.0f    %10.0f %10.0f %6.2f %10.1f\n", 
	   stat[ib].fBranch->GetName(),
	   stat[ib].XMean(), sqrt(stat[ib].SigXX()),
	   stat[ib].TotBytes(),
	   stat[ib].ZipBytes(),
	   stat[ib].CompRatio(),
	   100.0*diskFrac);
  }
  printf("-------------------------------------------------");
  printf("--------------------------------------------------\n");
  printf("........... total .............. <event size>");
  printf(" <sigma size> <disk size>\n"); 
  printf("-------------------------------------------------");
  printf("--------------------------------------------------\n");
  printf("%-30s %12.3f %12.3f %12.3f\n", " total event",
	 stat[nbranches].XMean(), sqrt(stat[nbranches].SigXX()),
	 meanCompSize);

  // Cleanup everything we booked here
  delete []stat;
  delete event;
}


//_____________________________________________________________________________
int TStnAna::AddArrays(TObjArray*  Arr1, TObjArray*  Arr2) {
  // in STNTUPLE files arrays are not allowed to have subfolders

//  TObjArrayIter it1(Arr1);
//  TObjArrayIter it2(Arr2);
//
//  TObject    *o1, *o2;
//  TH1        *h1, *h2; 
//  TObjArray  *a1, *a2;

//    while ((o1 = it1.Next())) {
//      o2 = Arr2->FindObject(o1->GetName());

//      if (o1->InheritsFrom("TH1")) {
//        h1   = (TH1*) o1;
//        h2   = (TH1*) o2;

//        *h1 = *h1 + *h2;
//      }
//      else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//  //-----------------------------------------------------------------------------
//  // lists
//  //-----------------------------------------------------------------------------
//        a1 = (TObjArray*) o1;
//        a2 = (TObjArray*) o2;
//        add_arrays(a1,a2);
//      }
//      else {
//        printf("%-30s %-20s *************** add_array EMOE *****************\n",
//  	     o1->GetName(),o1->ClassName());
//      }
//    }
  return 0;
}


//_____________________________________________________________________________
Int_t TStnAna::AddFolders(TFolder*   Fol1, TFolder*   Fol2) {

  // for each histogram from Fol1 add corresponding histogram from Fol2 to it

//  TIter it1(Fol1->GetListOfFolders());
//
//  TObject  *o1, *o2;
//  TH1      *h1, *h2; 
//  TFolder  *f1, *f2;
//  TObjArray *a1, *a2;
//-----------------------------------------------------------------------------
//  loop over the objects stored in this Fol1
//-----------------------------------------------------------------------------
//    while ((o1 = it1.Next())) {
//      o2 = Fol2->FindObject(o1->GetName());

//      if (o1->InheritsFrom("TH1")) {
//        h1 = (TH1*) o1;
//        h2 = (TH1*) o2;

//        *h1 = *h1 + *h2;
//      }
//      else if (strcmp(o1->ClassName(),"TFolder") == 0) {
//        printf("%-30s %-20s folder \n",o1->GetName(),o1->ClassName());
//        f1 = (TFolder*) o1;
//        f2 = (TFolder*) o2;
//  //-----------------------------------------------------------------------------
//  // compare histogram contents of 2 folders
//  //-----------------------------------------------------------------------------
//        AddFolders(f1,f2);
//      }
//      else if (strcmp(o1->ClassName(),"TObjArray") == 0) {
//  //-----------------------------------------------------------------------------
//  // compare lists of histograms
//  //-----------------------------------------------------------------------------
//        a1 = (TObjArray*) o1;
//        a2 = (TObjArray*) o2;
//        AddArrays(a1,a2);
//      }
//      else {
//        printf("%-30s %-20s *************** EMOE *****************\n",
//  	     o1->GetName(),o1->ClassName());
//      }
//    }
  return 0;
}

//_____________________________________________________________________________
int TStnAna::MergeHistograms(const char* List, const char* OutputFile) {
  // given List of STNTUPLE histogram files (i.e. "a/*.root") merges them and 
  // writes output into OutputFile


  char  cmd[200];

  //TFile*          f;
  TFile*          output_file;
  TFolder         *fol1, *fol2;

  int first = 1;

  char fn[200];
  FILE* file = gSystem->OpenPipe(Form("ls %s",List),"r");
  while ( fscanf(file,"%s",fn) > 0) {
    printf("%s\n",fn);
    if (first) {
      first = 0;
//-----------------------------------------------------------------------------
// figure out list of histograms
//-----------------------------------------------------------------------------
      sprintf(cmd,"cp %s %s",fn,OutputFile);
      gSystem->Exec(cmd);
      output_file = new TFile(OutputFile,"u");
    }
    else {
//-----------------------------------------------------------------------------
// loop over the histograms in the OutputFile and for each histogram add the
// same histogram from "fn" to it
//-----------------------------------------------------------------------------
//        fol1 = read_folder(OutputFile,"Ana");
//        fol2 = read_folder(fn,"Ana");
      AddFolders(fol1,fol2);
    }
  }
  fclose(file);
  output_file->Write();
  output_file->Close();

  return 0;
}

//_____________________________________________________________________________
void TStnAna::Clear(const char* Opt) {
  // assume TStnAna has been already initialized once, reinitialize
  // leave fMode as it was

  fInitialized      = 0;
				// delete contents of the lists, but not the
				// lists themselves
  fListOfRuns->Delete();
  fModuleList->Delete();
//-----------------------------------------------------------------------------
// so far deleting a folder leads to deletion of the histograms as well...
// ... hmm ... is Delete()
//-----------------------------------------------------------------------------
  fFolder->Clear();

  fEventNumber      = -1;
  fRunNumber        = -1;
  fMinRunNumber     = -1;
  fMaxRunNumber     = INT_MAX;
  fEvent            = 0;
  fEntry            = -1;
  fNProcessedEvents = 0;
  fNPassedEvents    = 0;
  fNEventsToReport  = 500;
  fOutputFile       = 0;
  fOutputTree       = 0;
  fOutputModule     = 0;
  fDBManager        = TStnDBManager::Instance();
  fGoodRun          = 1;
  fGoodRunList      = 0;
  fPrintLevel       = 0;

  TH1::AddDirectory(0);
}


//_____________________________________________________________________________
void TStnAna::Print(const char* Opt) const {
}
