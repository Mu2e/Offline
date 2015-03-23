//-----------------------------------------------------------------------------
// Description:
// -----------
// Class StntupleModule : base class for STNTUPLE modules
// It inherits from RootHistModule and adds static pointer to TStnEvent
// plus `AddDataBlock' method
//
// Nov 23 2000 P.Murat
//-----------------------------------------------------------------------------
#include <string>
#include <cstdio>

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TTree.h"
#include "TFolder.h"
#include "TROOT.h"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnErrorLogger.hh"
#include "Stntuple/obj/TStnDataBlock.hh"

#include "Stntuple/mod/StntupleModule.hh"

// ClassImp(StntupleModule)

TStnEvent*       StntupleModule::fgEvent          = 0;
TStnErrorLogger* StntupleModule::fgErrorLogger    = 0;
TFolder*         StntupleModule::fgStntupleFolder = 0;
//-----------------------------------------------------------------------------
// constructors
//-----------------------------------------------------------------------------
StntupleModule::StntupleModule(fhicl::ParameterSet const& PSet, const char* Name): 
  THistModule(PSet,Name) 
{
  if (! fgEvent      ) {
    fgEvent       = new TStnEvent();
    //    fgErrorLogger = new TStnErrorLogger();
    //    fgEvent->SetErrorLogger(fgErrorLogger);
    fgStntupleFolder = gROOT->GetRootFolder()->AddFolder("Stntuple",
							 "STNTUPLE folder");
    //    fgStntupleFolder->Add(fgErrorLogger);
    THistModule::fgMaxFileSize = 8000;
  }
}


//-----------------------------------------------------------------------------
StntupleModule::~StntupleModule() {
  // folders do not have to be deleted!
  if (fgEvent) {
    delete fgEvent;
    fgEvent = 0;
    //    delete fgErrorLogger;
    fgErrorLogger = 0;
  }
}

//_____________________________________________________________________________
void StntupleModule::LogError(const char* Message) 
{
  // remember that Message has format "Err<%i>: %s", so it is possible to 
  // retrieve the actual error code from it

  //  ERRLOG(ELwarning,Message) << endmsg;
}

//_____________________________________________________________________________
void StntupleModule::LogError(char* Message) 
{
  // remember that Message has format "Err<%i>: %s", so it is possible to 
  // retrieve the actual error code from it

  //  ERRLOG(ELwarning,Message) << endmsg;
  printf(">>> StntupleModule::LogError(char* Message): message: %s\n",Message);
}

//_____________________________________________________________________________
TStnDataBlock* 
StntupleModule::AddDataBlock(const char* branch_name,
			     const char* class_name,
			     Int_t       (*f)(TStnDataBlock*,AbsEvent*,Int_t),
			     Int_t       buffer_size,
			     Int_t       split,
			     Int_t       compression) 
{
  // adds new branch to fgTree and registers a data block corresponding to it

  TBranch*       branch;
  TStnDataBlock* block;

  int            rc;

  TStnNode*      node;
  node = 0;
  rc  = fgEvent->AddDataBlock(branch_name,class_name,node);

  if (rc == 0) {
				// everything is OK

    branch = fgTree->Branch(branch_name,class_name,
			    node->GetDataBlockAddress(),
			    buffer_size,
			    split);
    branch->SetCompressionLevel(compression);
    block = node->GetDataBlock();
    block->SetExternalInit(f);
    block->SetNode(node);
  }
  else if (rc > 0) {
				// already existing branch

    printf(" StntupleModule::AddDataBlock: ");
    printf(" an attempt to redefine the existing block made for");
    printf(" branch %s and class %s\n",branch_name,class_name);
    block = node->GetDataBlock();
  }
  else {
				// can't create branch

    printf(" StntupleModule::AddDataBlock : can\'t add block for");
    printf(" branch %s and class %s\n",branch_name,class_name);
    block = NULL;
  }
  return block;
}


//_____________________________________________________________________________
Int_t 
StntupleModule::SetResolveLinksMethod(const char* BlockName,
				      Int_t  (*f)(TStnDataBlock*,AbsEvent*,Int_t))
{
  // set ResolveLinks method for a data block BlockName

  TStnDataBlock* block;

  int rc = 0;

  block  = fgEvent->GetDataBlock(BlockName);

  if (block == 0) {
    rc = -1;
  }
  else {
    block->SetResolveLinksMethod(f);
  }
  return rc;
}








