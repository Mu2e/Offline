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
#include "TClass.h"

#include "Stntuple/obj/TStnErrorLogger.hh"

ClassImp(TStnErrorLogger)

//_____________________________________________________________________________
TStnErrorLogger::TStnErrorLogger(const char* Name, const char* Title):
  TNamed(Name,Title) {
}

//_____________________________________________________________________________
TStnErrorLogger::~TStnErrorLogger() {
}

//_____________________________________________________________________________
void TStnErrorLogger::Clear(Option_t* opt) {
}

//------------------------------------------------------------------------------
void TStnErrorLogger::Print(Option_t* opt) const {
  printf(" *** TStnErrorLogger::Print is not implemented yet\n");
}


//_____________________________________________________________________________
void TStnErrorLogger::Report(Int_t ErrorCode, const char* Message) {
  Emit("Report(Int_t,const char*)",Form("Err:%i: %s",ErrorCode,Message));
}
