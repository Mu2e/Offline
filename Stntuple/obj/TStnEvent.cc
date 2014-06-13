//-----------------------------------------------------------------------------
// this class is a run-time interface to the event data
// it keeps a list of the branches, requested by the various analysis modules 
// and is not supposed to be directly used by the unexperienced users
//-----------------------------------------------------------------------------

#ifdef __GNUG__
#pragma implementation
#endif
#include <math.h>
#include <iostream>
#include <iomanip>

#include "TROOT.h"
#include "TBranch.h"
#include "TClass.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnEvent.hh"

ClassImp(TStnEvent)
//_____________________________________________________________________________
TStnEvent::TStnEvent() {
  // list of nodes is the one which owns all of them
  
  fListOfNodes       = new TObjArray(10);
  fListOfInputNodes  = new TObjArray(10);
  fListOfOutputNodes = new TObjArray(10);
  fDropList          = new TObjArray(10);

  fListOfObjects     = new TObjArray(10);
  fListOfHptl        = new TObjArray(10);
  fListOfHptl->SetName("ListOfHptl");

  fListOfUnusedNodes = new TObjArray(10);

  fLastNumber = -1;
}


//_____________________________________________________________________________
TStnEvent::~TStnEvent()
{
//-----------------------------------------------------------------------------
//  input nodes are owned by the event itself
//-----------------------------------------------------------------------------
  //printf(" TStnEvent: input nodes.\n"); fflush(stdout); fflush(stderr);
  fListOfInputNodes->Delete();
  delete fListOfInputNodes;
//-----------------------------------------------------------------------------
//  nodes are owned by the event itself
//-----------------------------------------------------------------------------
  //printf(" TStnEvent: nodes.\n"); fflush(stdout); fflush(stderr);
  fListOfNodes->Clear();             // These are a subset of fListOfInputNodes
  delete fListOfNodes;
//-----------------------------------------------------------------------------
//  output nodes are owned by the event itself
//-----------------------------------------------------------------------------
  //printf(" TStnEvent: output nodes.\n"); fflush(stdout); fflush(stderr);
  fListOfOutputNodes->Clear();            // These are a subset of fListOfNodes
  delete fListOfOutputNodes;
//-----------------------------------------------------------------------------
//  list of unused nodes has to be deleted (not used as far as I see!)
//-----------------------------------------------------------------------------
  //printf(" TStnEvent: unused nodes.\n"); fflush(stdout); fflush(stderr);
  fListOfUnusedNodes->Delete();
  delete fListOfUnusedNodes;
//-----------------------------------------------------------------------------
//  all objects have their respective owners
//-----------------------------------------------------------------------------
  //printf(" TStnEvent: objects.\n"); fflush(stdout); fflush(stderr);
  fListOfObjects->Delete();
  delete fListOfObjects;

  //printf(" TStnEvent: hptl.\n"); fflush(stdout); fflush(stderr);
  fListOfHptl->Delete();
  delete fListOfHptl;

  //printf(" TStnEvent: drops.\n"); fflush(stdout); fflush(stderr);
  fDropList->Clear();
  delete fDropList;
}


//_____________________________________________________________________________
void TStnEvent::Clear(Option_t* opt) {
  // clear all the variables
  TStnNode* node;
  TIter     it(fListOfNodes);
  while ((node = ((TStnNode*) it.Next()))) {
    node->GetDataBlock()->Clear();
  }

  fListOfObjects->Clear();
  fListOfHptl->Clear();
}


//_____________________________________________________________________________
Int_t TStnEvent::ReadTreeEntry(Int_t Entry) {
  // read Entry from the current tree - kludge needed by the output module...
  // revisit... Assume that fCurrentEntry has already been set...
  int nb = 0;
  TIter     it(fListOfNodes);
  while (TStnNode* node = (TStnNode*) it.Next()) {
    nb += node->GetEntry(Entry);
  }
  return nb;
}

//_____________________________________________________________________________
Int_t TStnEvent::AddDataBlock(const char* branch_name, 
			      const char* class_name,
			      TStnNode*&  node) 
{
  // create new node for branch `branch_name' and C++ class `class_name'
  // class `class_name' should inherit from TStnDataBlock and have a 
  // ROOT dictionary with ClassDef/ClassImp

  TClass* cl = gROOT->GetClass(class_name);
  if (!cl) {
    Error("AddDataBlock","class %s doesn\'t have a dictionary",class_name);
    node = NULL;
    return -1;
  }

  if (! cl->InheritsFrom("TStnDataBlock")) {
    Error("AddDataBlock","class %s doesn\'t inherit from TStnDataBlock",
	  class_name);
    node = NULL;
    return -1;
  }
					// make sure there is no branch with 
					// the same name 

  node = (TStnNode*)  fListOfNodes->FindObject(branch_name);
  if (node) {
					// node already exists, generate 
					// warning
    return 1;
  }
					// everything is OK, create new node

  node = new TStnNode(branch_name,cl,this);
  fListOfNodes->Add(node);

  return 0;
}

//_____________________________________________________________________________
Int_t TStnEvent::AddOutputBlock(const char* BranchName, TStnDataBlock* Block) {
  // block already exists and its initialization is in the hands of user!

  TStnNode* node;
  Int_t    rc = 0;
					// make sure there is no branch with 
					// the same name 
  node = (TStnNode*)  fListOfOutputNodes->FindObject(BranchName);
  if (node) {
					// node already exists, report an error
    Error("AddDataBlock","block %s already registered",BranchName);
    return -1;
  }
					// everything is OK, create new node
					// and add it to the output list but 
					// not to the list of active input 
					// nodes
  node = new TStnNode(BranchName,Block,this);
  fListOfOutputNodes->Add(node);

  return rc;
}


//_____________________________________________________________________________
TStnDataBlock* TStnEvent::GetDataBlock(const char* branch_name) {
  // return pointer to the data block corresponding to branch `branch_name'

  TStnNode* node = (TStnNode*) fListOfNodes->FindObject(branch_name);
  				// try to find the block in the list of output
				// nodes 

  if (!node)
    node = (TStnNode*) fListOfOutputNodes->FindObject(branch_name);

  if (!node)
    return NULL;
  else
    return node->GetDataBlock();
}

//_____________________________________________________________________________
TStnDataBlock* TStnEvent::UnpackDataBlock(const char* branch_name) {

  TStnDataBlock* blk = GetDataBlock(branch_name);
  if(blk==NULL) return blk;

  int ientry = GetCurrentTreeEntry();
  blk->GetEntry(ientry);

  return blk;
}


//_____________________________________________________________________________
TStnDataBlock** TStnEvent::GetDataBlockAddress(const char* branch_name) {
  TStnNode* node = (TStnNode*) fListOfNodes->FindObject(branch_name);
  
  if (!node)
    node = (TStnNode*) fListOfOutputNodes->FindObject(branch_name);

  if (!node)
    return NULL;
  else
    return node->GetDataBlockAddress();
}


//_____________________________________________________________________________
void TStnEvent::SetEventNumber(Int_t RunNumber, 
			       Int_t EventNumber, 
			       Int_t SectionNumber) 
{
  fRunNumber     = RunNumber;
  fEventNumber   = EventNumber;
  fSectionNumber = SectionNumber;
}

//_____________________________________________________________________________
int TStnEvent::InitFromGenp() {
  // loop over the GENP particles and fill the lists of high-Pt electrons and
  // muons (no isolation  requirements for the moment, just Pt and eta cuts)
  return 0;
}

//______________________________________________________________________________
int TStnEvent::Init(AbsEvent* event, int mode) {

  Clear();

  TIter it(fListOfNodes);
  while (TStnNode* node = (TStnNode*) it.Next()) {
    TStnDataBlock* block = node->GetDataBlock();
//-----------------------------------------------------------------------------
// check if user initialization has been requested
//-----------------------------------------------------------------------------
    if (! block->UserInitialization()) {
      //printf("Calling init on block %s\n",node->GetName());
      block->Init(event,0);
    }
  }

  return 0;
}

//------------------------------------------------------------------------------
void TStnEvent::Print(Option_t* opt) const {
  TIter it(fListOfNodes);
  while (TStnNode* node = (TStnNode*) it.Next()) {
    printf(" ===================== data block %s",node->GetName());
    printf(" =====================\n");
    node->GetDataBlock()->Print(opt);
  }
}
