#include <iostream>
#include <iomanip>

#include "obj/TStnPidBlock.hh"
#include "obj/TStnPid.hh"

ClassImp(TStnPidBlock)

//______________________________________________________________________________
void TStnPidBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TStnPidBlock.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    R__b >> fNTracks;
    fListOfTrackPid->Streamer(R__b);
  } 
  else {
    R__b.WriteVersion(TStnPidBlock::IsA());
    R__b << fNTracks;
    fListOfTrackPid->Streamer(R__b);
  }
}

//_____________________________________________________________________________
TStnPidBlock::TStnPidBlock() {
  fNTracks   = 0;
  fListOfTrackPid = new TClonesArray("TStnPid",100);
  fListOfTrackPid->BypassStreamer(kFALSE);
  fCollName     = "default";
}


//_____________________________________________________________________________
TStnPidBlock::~TStnPidBlock() {
  fListOfTrackPid->Delete();
  delete fListOfTrackPid;
}


//_____________________________________________________________________________
void TStnPidBlock::Clear(Option_t* opt) {
  fNTracks = 0;
  fListOfTrackPid->Clear(opt);
}

//------------------------------------------------------------------------------
void TStnPidBlock::Print(Option_t* opt) const {

  int banner_printed = 0;
  for (int i=0; i<fNTracks; i++) {
    TStnPid* t = ((TStnPidBlock*) this)->Pid(i);
    if (! banner_printed) {
      t->Print("banner");
      banner_printed = 1;
    }
    t->Print("data");
  }
}

