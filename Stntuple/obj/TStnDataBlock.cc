//_____________________________________________________________________________
//  base class for STNTUPLE data block
// 
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 10 2000
//
//  TStnDataBlock's are responsible for memory cleanup, the corresponding
//  branches are created with fAutoDelete = 0;
//  by default data block is filled with the event data in BookStntupleModule,
//  if user initialization has been requested (by setting 
//  fInitializedByUser to 1),
//  BookStntupleModule doesn't fill this block. In this case it is a user
//  responsibility to fill the block
//
//  initialization of the data block with the event data is done either using
//  overloaded protected `fInit' function, or
//_____________________________________________________________________________
#include <iostream>

#include "TClass.h"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnEvent.hh"

const Float_t TStnDataBlock::kUndefined = 1.e6;

ClassImp(TStnDataBlock)
//_____________________________________________________________________________
TStnDataBlock::TStnDataBlock():
  fCollName("none")
{
  // by default initialization is done by the overloaded function

  fUserInitialization = 0;
  fExternalInit       = 0;
  fResolveLinks       = 0;
  fLinksInitialized   = 0;
  fInitMode           = 0;
  fNode               = 0;
  fEvent              = 0;
  fMessageList        = new TObjArray(10);
  fCurrentEntry       = -1;
  fValid              = 0;
  f_EventNumber       = -1;
  f_RunNumber         = -1;
  fListOfCollNames    = new TObjArray();
}


//_____________________________________________________________________________
TStnDataBlock::~TStnDataBlock() {
  fMessageList->Delete();
  delete fMessageList;
  delete fListOfCollNames;
}

//-----------------------------------------------------------------------------
// label coding: ModuleLabel@Description;Processname
//-----------------------------------------------------------------------------
TNamed*  TStnDataBlock::AddCollName(const char* CollName, 
				    const char* ModuleLabel, 
				    const char* Description, 
				    const char* ProcessName) {
  TNamed* n;

  n = (TNamed*) fListOfCollNames->FindObject(CollName);

  if (n == 0) {
    TString s = ModuleLabel;
    s        += '@';
    s        += Description;
    s        += ';';
    if (ProcessName[0] != 0) {
      s      += ProcessName;
    }
    n = new TNamed(CollName,s.Data());
    fListOfCollNames->Add(n);
  }
  
  return n;
}


//_____________________________________________________________________________
Int_t TStnDataBlock::GetEntry(Int_t TreeEntry) {
  // TreeEntry refers to the current entry in the tree, fEntry - in the chain

  int current_entry = fEvent->GetCurrentEntry();
  if (fCurrentEntry == current_entry) return 1;
  fCurrentEntry = current_entry;

  return fNode->GetEntry(TreeEntry);
}

//_____________________________________________________________________________
void TStnDataBlock::SetCollName(const char* Process, 
				const char* Description, 
				const char* CollType) {
  fCollName = Process;
  fCollName += "@";
  fCollName += Description;
  if (CollType) {
    fCollName += ";";
    fCollName += CollType;
  }
}

//-----------------------------------------------------------------------------
void TStnDataBlock::GetModuleLabel(const char* CollName, char* ModuleLabel) {

  int len(0);

  TNamed*  n = (TNamed*) fListOfCollNames->FindObject(CollName);

  if (n != NULL) {
    TString s(n->GetTitle());

    len = s.Index('@');
    strncpy(ModuleLabel,s.Data(),len);
    ModuleLabel[len] = 0;

    //    printf(" collname = %-40s  n = %-s\n",CollName,n->GetTitle());
  }
  else {
    ModuleLabel[0] = 0;
    //    printf(" collname = %-40s  n = ... undefined...\n",CollName);
  }
}

//_____________________________________________________________________________
void TStnDataBlock::GetDescription(const char* CollName, char* Description) {
  int len(0);
  TNamed*  n = (TNamed*) fListOfCollNames->FindObject(CollName);

  if (n) {
    TString s(n->GetTitle());

    int first = s.Index('@')+1;
    len = s.Index(';')-first;
    strncpy(Description,s.Data()+first,len);
  }
  Description[len] = 0;
}

//-----------------------------------------------------------------------------
// string of a form: MMMMM@DDDDDDD;PPPPPP
//
// MMMMMM : module name
// DDDDDD : description
// PPPPPP : process name
//-----------------------------------------------------------------------------
void TStnDataBlock::GetProcessName(const char* CollName, char* ProcessName) {

  int first, nc, len(0);

  TNamed*  n = (TNamed*) fListOfCollNames->FindObject(CollName);

  if (n != NULL) {
    TString s(n->GetTitle());

    first = s.Index(';')+1;
    len   = s.Length();
    nc    = len-first;
    if (nc > 0) strncpy(ProcessName,s.Data()+first,nc);
    ProcessName[nc] = 0;
  }
  else {
    ProcessName[0] = 0;
  }
}

//_____________________________________________________________________________
void TStnDataBlock::Print(Option_t* option) const {
  std::cout << "--------------------------------------" << std::endl;
  std::cout << "TStnDataBlock Dump for " << fCollName <<std::endl;
  std::cout << "Run="<<f_RunNumber<<"  Event="<<f_EventNumber<<std::endl;
  std::cout << "Links init="<<fLinksInitialized 
	    << "   User init="<< fUserInitialization<<std::endl;
  std::cout << "Current Entry="<<fCurrentEntry 
	    << "Valid" << fValid<<std::endl;

}




