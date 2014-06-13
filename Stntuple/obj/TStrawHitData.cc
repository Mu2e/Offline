///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TStrawHitData
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/obj/TStrawHitData.hh"

ClassImp(TStrawHitData)

//_____________________________________________________________________________
void TStrawHitData::Streamer(TBuffer &R__b) {
  int nwf = 3;
  
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    R__b >> fIndex;
    R__b.ReadFastArray(&fTime,nwf);
  }
  else {
    R__b.WriteVersion(TStrawHitData::IsA());
    R__b << fIndex;
    R__b.WriteFastArray(&fTime,nwf);
  } 
}

//_____________________________________________________________________________
TStrawHitData::TStrawHitData(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TStrawHitData::~TStrawHitData() {
}

//_____________________________________________________________________________
void TStrawHitData::Clear(Option_t* opt) {
  fIndex  = -1;
  fTime   = 1.e6;
  fDt     = 1.e6;
  fEnergy = -1;
}

//_____________________________________________________________________________
void TStrawHitData::Print(Option_t* opt) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
}
