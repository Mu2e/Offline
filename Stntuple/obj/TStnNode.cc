//-----------------------------------------------------------------------------
//  base class for STNTUPLE node (holder of TStnDataBlock)
// 
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 10 2000
//-----------------------------------------------------------------------------
#include "TBranch.h"
#include "TClass.h"

#include <Stntuple/obj/TStnNode.hh>
#include <Stntuple/obj/TStnDataBlock.hh>

ClassImp(TStnNode)

//_____________________________________________________________________________
TStnNode::TStnNode() {
  fObject       = 0;
  fDeleteObject = 0;
}


//_____________________________________________________________________________
TStnNode::TStnNode(const char* BranchName, 
		   TClass*     Class,
		   TStnEvent*  Event,
		   Int_t       (*F)(TStnDataBlock *,TStnEvent *,Int_t)):
  TNamed(BranchName,BranchName)
{
  // one needs to make sure that cl inherits from TStnDataBlock before 
  // creating a node
  fBranch       = NULL;
  fEvent        = Event;
  fObject       = (TStnDataBlock*) Class->New();
  fObject->SetEvent(Event);
  fFunc         = F;
  fDeleteObject = 1;
}


//_____________________________________________________________________________
TStnNode::TStnNode(TBranch*    Branch, 
		   TClass*     Class,
		   TStnEvent*  Event,
		   Int_t       (*F)(TStnDataBlock *,TStnEvent *,Int_t)):
  TNamed(Branch->GetName(),Branch->GetName())
{
  // one needs to make sure that cl inherits from TStnDataBlock before 
  // creating a node
  fBranch       = Branch;
  fEvent        = Event;
  fObject       = (TStnDataBlock*) Class->New();
  fObject->SetEvent(Event);
  fFunc         = F;
  fDeleteObject = 1;
}


//_____________________________________________________________________________
TStnNode::TStnNode(const char*    BranchName, 
		   TStnDataBlock* Block,
		   TStnEvent*     Event,
		   Int_t          (*F)(TStnDataBlock *,TStnEvent *,Int_t)):
  TNamed(BranchName,BranchName)
{
  // one needs to make sure that cl inherits from TStnDataBlock before 
  // creating a node

  fBranch = 0;
  fEvent  = Event;
  fObject = Block;
  fObject->SetEvent(Event);
  fFunc   = F;
}


//_____________________________________________________________________________
TStnNode::~TStnNode() {
  if (fDeleteObject)
    delete fObject;
}


//_____________________________________________________________________________
Int_t TStnNode::GetEntry(Int_t Ientry) {
  // fBranch = 0 means that this node is unused

  if (fFunc) {
    return fFunc(fObject,fEvent,0);
  }
  else if (fBranch) {
    return fBranch->GetEntry(Ientry);
  }
  else {
    return -1;
  }
}

//_____________________________________________________________________________
void TStnNode::Print(Option_t* option) const {
}
