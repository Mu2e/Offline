//
// 2014-07-15 G. Pezzullo
//

#include "Stntuple/obj/TVdetDataBlock.hh"

ClassImp(TVdetDataBlock)

//-----------------------------------------------------------------------------
  void TVdetDataBlock::Streamer(TBuffer &R__b) {
  if(R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
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
void TVdetDataBlock::Print(Option_t* Option) const {
  // print all hits in the virtual detectors

  int banner_printed(0);
  const TVdetHitData* hit;

  TString opt = Option;
  opt.ToLower();

  for(int i=0; i<fNHits; i++) {
    hit = ((TVdetDataBlock*) this)->Hit(i);
    if ((opt == "") || (opt.Index("banner") >= 0)) {
      if (banner_printed == 0) {
	hit->Print("banner");
	banner_printed = 1;
      }
    }

    if ((opt == "") || (opt.Index("data") >= 0)) {
      hit->Print("data");
    }
  }
}

