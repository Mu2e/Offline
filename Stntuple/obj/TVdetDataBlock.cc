//
// 2014-07-15 G. Pezzullo
//

#include "Stntuple/obj/TVdetDataBlock.hh"

ClassImp(TVdetDataBlock)

//-----------------------------------------------------------------------------
  void TVdetDataBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    R__b >> fNHits;
    fListOfHits->Streamer(R__b);
  }
  else {
    R__b.WriteVersion(TVdetDataBlock::IsA());
    R__b << fNHits;
    fListOfHits->Streamer(R__b);
  }
}

//______________________________________________________________________________
TVdetDataBlock::TVdetDataBlock() {
  fListOfHits = new TClonesArray("TVdetHitData",30240);
  fListOfHits->BypassStreamer(kFALSE);
  Clear();
}

//______________________________________________________________________________
TVdetDataBlock::~TVdetDataBlock() {
  fListOfHits->Delete();
  delete fListOfHits;
}

//______________________________________________________________________________
void TVdetDataBlock::Clear(Option_t* opt) {
  fListOfHits->Clear();
  fNHits=0;
}

//______________________________________________________________________________
void TVdetDataBlock::Print(Option_t* opt) const {
  // print all hits in the virtual detectors
  printf(" *** print Virtual detectors hits *** \nNumber of hits: %d\n",fNHits);
  for(int i=0; i<fNHits; i++) {
    fListOfHits->At(i)->Print();
  }
}

