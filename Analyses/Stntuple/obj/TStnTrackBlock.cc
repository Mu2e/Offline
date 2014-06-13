#include <iostream>
#include <iomanip>

#include "obj/TStnTrackBlock.hh"
#include "obj/TStnTrack.hh"

ClassImp(TStnTrackBlock)

//______________________________________________________________________________
void TStnTrackBlock::Streamer(TBuffer &R__b) {
  // Stream an object of class TStnTrackBlock.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    R__b >> fNTracks;
    fListOfTracks->Streamer(R__b);
    for (int i=0; i<fNTracks; i++) {
      Track(i)->SetNumber(i);
    }
  } 
  else {
    R__b.WriteVersion(TStnTrackBlock::IsA());
    R__b << fNTracks;
    fListOfTracks->Streamer(R__b);
  }
}

//_____________________________________________________________________________
TStnTrackBlock::TStnTrackBlock() {
  fNTracks   = 0;
  fListOfTracks = new TClonesArray("TStnTrack",100);
  fListOfTracks->BypassStreamer(kFALSE);
  fCollName     = "default";
}


//_____________________________________________________________________________
TStnTrackBlock::~TStnTrackBlock() {
  fListOfTracks->Delete();
  delete fListOfTracks;
}


//_____________________________________________________________________________
void TStnTrackBlock::Clear(Option_t* opt) {
  fNTracks = 0;
  fListOfTracks->Clear(opt);
}

//------------------------------------------------------------------------------
void TStnTrackBlock::Print(Option_t* opt) const {

  int banner_printed = 0;
  for (int i=0; i<fNTracks; i++) {
    TStnTrack* t = ((TStnTrackBlock*) this)->Track(i);
    if (! banner_printed) {
      t->Print("banner");
      banner_printed = 1;
    }
    t->Print("data");
  }
}

