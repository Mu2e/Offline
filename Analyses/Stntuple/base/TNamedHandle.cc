///////////////////////////////////////////////////////////////////////////////
// handle doesn't own the object.
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/base/TNamedHandle.hh"


ClassImp(TNamedHandle)

//-----------------------------------------------------------------------------
TNamedHandle::TNamedHandle(): TNamed() {
  fObject = NULL;
}


//-----------------------------------------------------------------------------
TNamedHandle::TNamedHandle(const char* Name, void* Object): TNamed(Name,Name) {
  fObject = Object;
}


//_____________________________________________________________________________
TNamedHandle::~TNamedHandle() {
}


// } // end namespace




