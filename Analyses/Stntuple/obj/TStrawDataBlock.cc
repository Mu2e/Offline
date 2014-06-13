//
// 2014-01-27 P.Murat
//

#include "Stntuple/obj/TStrawDataBlock.hh"

ClassImp(TStrawDataBlock)

//-----------------------------------------------------------------------------
  void TStrawDataBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    R__b >> fNHits;
    fListOfHits->Streamer(R__b);
  }
  else {
    R__b.WriteVersion(TStrawDataBlock::IsA());
    R__b << fNHits;
    fListOfHits->Streamer(R__b);
  }
}

//______________________________________________________________________________
TStrawDataBlock::TStrawDataBlock() {
  fListOfHits = new TClonesArray("TStrawHitData",30240);
  fListOfHits->BypassStreamer(kFALSE);
  Clear();
}

//______________________________________________________________________________
TStrawDataBlock::~TStrawDataBlock() {
  fListOfHits->Delete();
  delete fListOfHits;
}

//______________________________________________________________________________
void TStrawDataBlock::Clear(Option_t* opt) {
  fListOfHits->Clear();
  fNHits=0;
}

//______________________________________________________________________________
void TStrawDataBlock::Print(Option_t* opt) const {
  // print all hits in the straw tracker
  printf(" *** print Straw tracker *** \nNumber of hits: %d\n",fNHits);
  for(int i=0; i<fNHits; i++) {
    fListOfHits->At(i)->Print();
  }
}

