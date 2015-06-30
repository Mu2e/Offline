//-----------------------------------------------------------------------------
//  Dec 28 2000 P.Murat: base class for STNTUPLE input module
//-----------------------------------------------------------------------------
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TString.h"

#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/obj/TStnEvent.hh"

ClassImp(TStnInputModule)

//_____________________________________________________________________________
TStnInputModule::TStnInputModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fDatasetList = new TList();
  fOwnChain    = false;
  fChain       = 0;
  fSplitInd = 0;
  fSplitTot = 0;
  fNEntries = 0;
  fNFiles   = 0;
}

//_____________________________________________________________________________
TStnInputModule::~TStnInputModule()
{
  // Get ride of the dataset list
  TListIter it(fDatasetList);
  while (TStnDataset* ds = (TStnDataset*) it.Next())
    delete ds;
  delete fDatasetList;
  // Sometimes the chain is created by us so be careful who owns it
  if (fOwnChain)
    delete fChain;
}

//_____________________________________________________________________________
int TStnInputModule::InitChain(const char* FileName, const char* TreeName)
{
  // assume that all the files in the directory should be chained
  void*   dir;
  TString filename(FileName);
 
  fOwnChain = true;
  fChain    = new TChain(TreeName);

  if (FileName == 0) {
    Warning("InitChain","input file is not defined -> Empty Chain");
    return -1;
  }
  // if this is not a remote file, and is a directory
  if ( TString(FileName).Index(":")<0 && (dir = gSystem->OpenDirectory(FileName))) {
    TString dirName = filename;
    const char* fn;
    while ((fn = gSystem->GetDirEntry(dir))) {
      if ((strcmp(fn,".") != 0) && (strcmp(fn,"..") != 0)) {
	filename = dirName + TString("/") + TString(fn);
	fChain->AddFile(filename.Data(),TChain::kBigNumber);
      }
    }
  }
  else {
				// this could be a single file
    // This is a single file ..
    fChain->AddFile(filename.Data(),TChain::kBigNumber);
  }
  if (fPrintLevel != 0)
    fChain->GetListOfFiles()->Print();

  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::AddDataset(TStnDataset* Dataset, int Print)
{

  fDatasetList->Add(Dataset);

  if (Print>0) {
    printf(" TStnInputModule::AddDataset - Adding \"%s\"\n",Dataset->GetName());
    fChain->GetListOfFiles()->Print();
  }

  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::BeginJob() {
  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::BeginRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::Event(Int_t i) {
  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::EndRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnInputModule::EndJob() {
  return 0;
}

//_____________________________________________________________________________
TStnNode* TStnInputModule::GetNode(const char* BranchName,  
				   const char* ClassName) {
  return 0;
}


//_____________________________________________________________________________
Int_t TStnInputModule::LoadEntry(Int_t Entry) {
  // deferred input: set the environment to read one entry

  if (! fChain)
    return -5;

  Int_t centry = fChain->LoadTree(Entry);
  if (centry < 0) {
    if (fPrintLevel > 0)
      printf(" TStnInputModule::LoadEntry - LoadTree = %d loading Entry: %d\n",
	     centry,Entry);
    return centry;
  }

  Int_t itree = fChain->GetTreeNumber();
  if (itree < 0 ) {
    printf("Error: TStnInputModule::LoadEntry - chain did not contain any trees\n");
    return -1;
  }
  if (itree != fCurrent) {
    fCurrent = itree;
    SetBranches();
    // Say we opened a new file and flush the output
    if (PrintLevel() == 1) {
      printf(" Opened new file: %s\n",fChain->GetFile()->GetName());
      fflush(stdout); fflush(stderr);
    }
  }

  GetEvent()->SetCurrentEntry(Entry);
  GetEvent()->SetCurrentTreeEntry(centry);

  return centry;
}

//_____________________________________________________________________________
int TStnInputModule::NextEvent(Int_t IEntry) {
  return LoadEntry(IEntry);
}

//_____________________________________________________________________________
int TStnInputModule::FindEvent(Int_t Run, Int_t Event) {
  return 0;
}
