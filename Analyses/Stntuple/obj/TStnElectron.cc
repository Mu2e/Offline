// P.Murat 
#include <iostream>

#include "obj/TStnElectron.hh"
#include "obj/TStnTrack.hh"
#include "obj/TStnCluster.hh"

ClassImp(TStnElectron)

//-----------------------------------------------------------------------------
void TStnElectron::Streamer(TBuffer& R__b) {
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
  }
  else {
    R__b.WriteVersion(TStnElectron::IsA());
  }
}

//_____________________________________________________________________________
TStnElectron::TStnElectron(Int_t Number) : TObject () {
    // 'Number' can be -1 ...

  fNumber = Number;
  
  fCluster = 0;
  fTrack   = 0;
}


//_____________________________________________________________________________
TStnElectron::~TStnElectron() {
}


//_____________________________________________________________________________
void TStnElectron::Print(Option_t* opt) const {
  printf("TStnElectron::Print <WARNING> Not implemented yet\n");
}
