//-----------------------------------------------------------------------------
//  Dec 28 2000 P.Murat: Run II STNTUPLE input module
//-----------------------------------------------------------------------------
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "TSystem.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TLeafObject.h"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TStnRun2InputModule.hh"

ClassImp(TStnRun2InputModule)
//_____________________________________________________________________________
TStnRun2InputModule::TStnRun2InputModule(): TStnInputModule() {
}

//_____________________________________________________________________________
TStnRun2InputModule::TStnRun2InputModule(const char* FileName, 
					 const char* TreeName):
  TStnInputModule("Run2InputModule","Run II Input Module")
{
  InitChain(FileName,TreeName);
  fCurrent = -1;
}


//_____________________________________________________________________________
TStnRun2InputModule::TStnRun2InputModule(TChain* chain):
  TStnInputModule("Run2InputModule","Run II Input Module")
{
  fChain   = chain;
  fCurrent = -1;
}


//_____________________________________________________________________________
TStnRun2InputModule::TStnRun2InputModule(TStnDataset* Dataset):
  TStnInputModule("Run2InputModule","Run II Input Module")
{
  fDatasetList->Add(Dataset);
  fChain   = Dataset->GetChain();
  fCurrent = -1;
}


//_____________________________________________________________________________
TStnRun2InputModule::~TStnRun2InputModule()
{
}

//_____________________________________________________________________________
Int_t TStnRun2InputModule::RegisterInputBranches(TStnEvent* Event) 
{
  // for each input branch create a node and add it to the ListOfNodes
  // this registration doesn't mean that all the branches will be read in 
  TObjArray* list = fChain->GetListOfBranches();
  TIter it(list);

  TBranch*          input_branch;
  TLeafObject*      leaf;
  const char*       class_name;
  TStnNode*         node;
  TClass*           cl;

  TObjArray* list_of_nodes = Event->GetListOfInputNodes();
  // Should be empty to begin with and analysis is responsible anyways
  list_of_nodes->Print();
  //list_of_nodes->Delete();

  while ((input_branch = (TBranchElement*) it.Next())) {
				// branch can be either TBranchElement or
				// TBranchObject
    if (strcmp(input_branch->ClassName(),"TBranchElement") == 0)
      class_name = ((TBranchElement*)input_branch)->GetClassName();
    else {
      leaf = (TLeafObject*)
	((TBranchObject*)input_branch)->GetLeaf(input_branch->GetName());
      class_name = leaf->GetTypeName();
    }

    cl = gROOT->GetClass(class_name);
    if (!cl || !cl->InheritsFrom("TStnDataBlock"))
      Error("RegisterInputBranches",
	    "class %s does not inheriting from TStnDataBlock",class_name);
    else {
      node = new TStnNode(input_branch,cl,Event);
      list_of_nodes->Add(node);
    }
  }
  return 0;
}


//_____________________________________________________________________________
TStnNode* TStnRun2InputModule::GetNode(const char* BranchName,
				       const char* ClassName)
{
				// check if a node with the name exists
  TStnNode* node;

  TStnEvent* ev = fAna->GetEvent();
  node = ev->FindNode(BranchName);

  if (node) return node;
				// node doesn't exist,
				// make sure branch exists

  TBranch* b = fChain->GetBranch(BranchName);
  if (!b) {
    Error("GetNode","branch %s doesn\'t exist",BranchName);
    return NULL;
  }

  TClass* cl = gROOT->GetClass(ClassName);
  if (!cl || !cl->InheritsFrom("TStnDataBlock")) {
    Error("GetNode","wrong class name %s",ClassName);
    return NULL;
  }
				// 

  node         = new TStnNode(BranchName,cl,ev);
  ev->GetListOfNodes()->Add(node);
  return node;
}

//_____________________________________________________________________________
Int_t TStnRun2InputModule::SetBranches() 
{
  // called by LoadTree when loading new file to get branch pointers

  TStnEvent* ev = fAna->GetEvent();

  TIter it(ev->GetListOfNodes());

  while (TStnNode* node = (TStnNode*) it.Next()) {
    const char* branch_name = node->GetName();
    TBranch* b = fChain->GetBranch(branch_name);
    node->SetBranch(b);
    if (b != 0) {
      b->SetAddress(node->GetDataBlockAddress());
      b->SetAutoDelete(0);
    }
    else {
      Error("SetBranches",Form("%s doesn\'t have branch %s",
			       GetChain()->GetFile()->GetName(),
			       branch_name));
    }
  }
  return 0;
}


//_____________________________________________________________________________
int TStnRun2InputModule::BeginJob() {

  fNEntries = 0;
  fNFiles = 0;

  if (fDatasetList->GetSize() > 0) {
    fOwnChain = true;

    TObjArray* arr = new TObjArray(0);
    arr->Clear();

    // loop over datasets
    int ntotev = 0;
    TListIter ids(fDatasetList);
    while (TStnDataset* ds = (TStnDataset*) ids.Next()) {
      TObjArrayIter it(ds->GetChain()->GetListOfFiles());
      // loop over files in each dataset
      //ds->GetChain()->GetListOfFiles()->Print();
      while (TChainElement* ce = (TChainElement*) it.Next()) {
	arr->AddLast(ce);
	ntotev += ce->GetEntries();
      }
    }
    int ntot = arr->GetEntries();

    fChain = new TChain("STNTUPLE");
    fCurrent = -1;

    // find pointers to the first file (counting from 0)
    // and the last plus 1, by splitting as requested
    if(ntot == 0 ) {
      printf("TStnRun2InputModule::BeginJob: Error: no files in dataset\n");
      return 1; 
    }

    // if no splitting, then run all files
    int ifirst = 0;
    int ilastp1 = ntot;

    if(fSplitTot>0) {
      printf("TStnRun2InputModule::BeginJob: datasets have %4d files %9d events\n",
	   ntot,ntotev);
      if(fSplitInd<0) {
	printf("TStnRun2InputModule::BeginJob: Error: SplitInd<0\n");
	return 1;
      }
      if(ntot<=fSplitTot) {
	ifirst = fSplitInd;
	ilastp1 = fSplitInd + 1;
	if(fSplitInd>=ntot) {
	  printf("TStnRun2InputModule::BeginJob: Error: no files for this job after split\n");
	  return 1;
	}
      } else { // ntot>fSplitTot
	if(fSplitInd>=fSplitTot) {
	  printf("TStnInputModule::BeginJob: Error: SplitInd>=SplitTot\n");
	  return 1;
	}
	int del = ntot/fSplitTot+1;
	int nbig = ntot - (del - 1)*fSplitTot;
	if(fSplitInd<nbig) {
	  ifirst = del*fSplitInd;
	  ilastp1 = ifirst + del;
	} else {
	  ifirst = del*nbig + (del-1)*(fSplitInd-nbig);
	  ilastp1 = ifirst + del - 1;
	}
      }
      printf("TStnRun2InputModule::BeginJob: Splitting, running on files %d to %d\n",ifirst,ilastp1-1);
    }

    fChain = new TChain("STNTUPLE");
    TObjArrayIter it(arr);
    int ind = 0;
    while (TChainElement* ce = (TChainElement*) it.Next()) {
      if(ind>=ifirst && ind<ilastp1) {
	fChain->AddFile(ce->GetTitle(),ce->GetEntries());
	fNEntries += ce->GetEntries();
	fNFiles++;
      }
      ind++;
    }
  } else {
    printf("TStnRun2InputModule::BeginJob Warning - no metadata,\n     opening all chained files to count entries...\n");
    fNEntries = int(fChain->GetEntries());
    fNFiles = fChain->GetListOfFiles()->GetEntries();
  }

  printf("TStnRun2InputModule::BeginJob: chained %4d files, %9d events\n",
	 fNFiles,fNEntries);

  fChain->UseCache(10);

  if (fAna->GetEvent() == 0) {
    fAna->SetEvent(new TStnEvent());
  }
  return 0;
}

//_____________________________________________________________________________
int TStnRun2InputModule::BeginRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnRun2InputModule::Event(Int_t i) {
  return 0;
}

//_____________________________________________________________________________
int TStnRun2InputModule::EndRun() {
  return 0;
}

//_____________________________________________________________________________
int TStnRun2InputModule::EndJob() {
  return 0;
}
