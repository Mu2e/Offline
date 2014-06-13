///////////////////////////////////////////////////////////////////////////////
//  Apr 28 2001 P.Murat: simple implementation for DB-type utilities
//  need to know MYRON flag for a given run anyway...
///////////////////////////////////////////////////////////////////////////////
#include "TDirectory.h"

// #include "Stntuple/data/TStnBeamPos.hh"
// #include "Stntuple/data/TStnBeamPosBlock.hh"
#include "Stntuple/obj/TStnRunSummary.hh"
// #include "Stntuple/obj/TStnTriggerTable.hh"
#include "Stntuple/obj/TStnDBManager.hh"
// #include "Stntuple/obj/TStnDeadList.hh"

#include <iostream>
#include <cstdio>
ClassImp(TStnDBManager)

TStnDBManager*  TStnDBManager::fgInstance = 0;

//_____________________________________________________________________________
TStnDBManager::TStnDBManager() 
{
  // Create the list of objects which we keep track of
  fListOfDbObjects = new TList();

  // Add some standard run dependent objects
  fListOfDbObjects->Add(new TStnRunSummary());
  // fListOfDbObjects->Add(new TStnTriggerTable());
  // fListOfDbObjects->Add(new TStnDeadList());

  // Beamline in the new scheme (time dependent)
  // fListOfDbObjects->Add(new TStnBeamPosBlock("CotBeamPosBlock"));
  // fListOfDbObjects->Add(new TStnBeamPosBlock("SvxBeamPosBlock"));

  //----------------------------------------------------------------------------
  // For backward compatibility
  //----------------------------------------------------------------------------
  // TStnBeamPos* bp;
  // bp = new TStnBeamPos("CotBeamPos");
  // bp->SetBit(kInvalidObject);      // Make sure this friend does not get written
  // fListOfDbObjects->Add(bp);
  // bp = new TStnBeamPos("SvxBeamPos");
  // bp->SetBit(kInvalidObject);      // Make sure this friend does not get written
  // fListOfDbObjects->Add(bp);
}

//_____________________________________________________________________________
TStnDBManager::~TStnDBManager() {
  fListOfDbObjects->Delete();
}

//_____________________________________________________________________________
TStnDBManager*  TStnDBManager::Instance() { 
  static Cleaner cleaner;
  return (fgInstance) ? fgInstance : (fgInstance = new TStnDBManager());
}

//------------------------------------------------------------------------------
TStnDBManager::Cleaner::Cleaner() {
}

//------------------------------------------------------------------------------
TStnDBManager::Cleaner::~Cleaner() {
  if (TStnDBManager::fgInstance) {
    delete TStnDBManager::fgInstance;
    TStnDBManager::fgInstance = 0;
  }
}


//_____________________________________________________________________________
TObject* TStnDBManager::GetTable(const char* Name)
{
  // Catch special case
  if (! strcmp("TStnBeamPos",Name)) 
    return fListOfDbObjects->FindObject("CotBeamPos");

  // Standard search
  TObject* obj = fListOfDbObjects->FindObject(Name);
  if (obj)
    return obj;

  // So far the object has not been found. If the object requested was
  // the shorthand name like "RunSummary" instead of "TStnRunSummary",
  // try to look for the latter
  TString TStnName(Name);
  TStnName.Prepend("TStn");
  obj = fListOfDbObjects->FindObject(TStnName.Data());
  if (! obj)
    printf(" Object  %s  does not exist!!!! -> CRASH ??\n",Name);

  return obj;
}


//_____________________________________________________________________________
Int_t TStnDBManager::Write(const char* Name, Int_t Option, Int_t BufSize) 
{
  int rc = 0;
  TIter it(fListOfDbObjects);
  while (TObject* obj = it.Next()) {
    if (! obj->TestBit(kInvalidObject))
      rc += obj->Write();
  }
  return rc;
}

//_____________________________________________________________________________
Int_t TStnDBManager::Read(const char* Name) 
{
  int nb = 0;
  TIter it(fListOfDbObjects);
  while (TObject* obj = it.Next()) {
    if (gDirectory->GetListOfKeys()->FindObject(obj->GetName()) != 0) {
      nb += obj->Read(obj->GetName());
    }
    else
      obj->Clear();
  }

  //----------------------------------------------------------------------------
  // Fill SvxBeamPos and CotBeamPos on the fly, if they are not available
  //----------------------------------------------------------------------------
  //  MakeAvgBeamline();

  //----------------------------------------------------------------------------
  // For compatibility with old file BeamPos: look for CotBeam only
  //----------------------------------------------------------------------------
  // if (gDirectory->GetListOfKeys()->FindObject("TStnBeamPos") != 0) {
  //   TStnBeamPos* cot = (TStnBeamPos*) GetTable("CotBeamPos");
  //   nb += cot->Read("TStnBeamPos");
  // }

  return nb;
}

