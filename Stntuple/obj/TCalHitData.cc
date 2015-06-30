///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TCalHitData
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/obj/TCalHitData.hh"

ClassImp(TCalHitData)

//_____________________________________________________________________________
void TCalHitData::Streamer(TBuffer &R__b) {
  int nwi = 2;
  int nwf = 2;
  
  if(R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion();
    R__b.ReadFastArray(&fID  ,nwi);
    R__b.ReadFastArray(&fTime,nwf);
  }
  else {
    R__b.WriteVersion(TCalHitData::IsA());
    R__b.WriteFastArray(&fID  ,nwi);
    R__b.WriteFastArray(&fTime,nwf);
  } 
}

//_____________________________________________________________________________
TCalHitData::TCalHitData(): TObject() {
  Clear();
}

//_____________________________________________________________________________
TCalHitData::~TCalHitData() {
}

//_____________________________________________________________________________
void TCalHitData::Clear(Option_t* opt) {
  fID        = -1;
  fNChannels = -1;
  fTime      = -1.;
  fEnergy    = -1.;
}

//_____________________________________________________________________________
void TCalHitData::Print(Option_t* opt) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
}
