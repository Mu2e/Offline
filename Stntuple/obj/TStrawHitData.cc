///////////////////////////////////////////////////////////////////////////////
//  2014-01-26 P.Murat TStrawHitData
///////////////////////////////////////////////////////////////////////////////
#include "TString.h"

#include "Stntuple/obj/TStrawHitData.hh"

ClassImp(TStrawHitData)

//-----------------------------------------------------------------------------
void TStrawHitData::ReadV1(TBuffer &R__b) {
  struct TStrawHitDataV1_t {
    int     fIndex;
    float   fTime;
    float   fDt;
    float   fEnergy;
  };

  TStrawHitDataV1_t data;

  int nwf_v1 = 3;

  R__b >> fIndex;
  R__b.ReadFastArray(&data.fTime,nwf_v1);

  fTime   = data.fTime;
  fDt     = data.fDt;
  fEnergy = data.fEnergy;
}



//_____________________________________________________________________________
void TStrawHitData::Streamer(TBuffer &R__b) {

  int nwi = ((int*) &fTime) - &fIndex;
  int nwf = &fMcMomentum - &fTime +1;
  
  if(R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion();
    if (R__v == 1) ReadV1(R__b);
    else {
//-----------------------------------------------------------------------------
// curent version: V2
//-----------------------------------------------------------------------------
      R__b.ReadFastArray(&fIndex,nwi);
      R__b.ReadFastArray(&fTime ,nwf);
    }
  }
  else {
    R__b.WriteVersion(TStrawHitData::IsA());
    R__b.WriteFastArray(&fIndex,nwi);
    R__b.WriteFastArray(&fTime ,nwf);
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
void TStrawHitData::Set(int Index, float Time, float Dt, float EnergyDep,
			int PdgCode, int MotherPdgCode, int GenCode, 
			int SimID, float McMomentum) 
{
  fIndex         = Index; 
  fTime          = Time; 
  fDt            = Dt; 
  fEnergy        = EnergyDep;
  fPdgCode       = PdgCode;
  fMotherPdgCode = MotherPdgCode;
  fGeneratorCode = GenCode;
  fSimID         = SimID;
  fMcMomentum    = McMomentum;
  
}

//_____________________________________________________________________________
void TStrawHitData::Clear(Option_t* opt) {
  fIndex         = -1;
  fTime          = 1.e6;
  fDt            = 1.e6;
  fEnergy        = -1;
  fPdgCode       = -1;
  fMotherPdgCode = -1;
  fGeneratorCode = -1;
  fSimID         = -1;
  fMcMomentum    = -1;
}

//_____________________________________________________________________________
void TStrawHitData::Print(Option_t* Option) const {
  // print Straw hit properties
  //  printf("Superlayer: %d, Wire: %d, Cell: %d,\n",fSuperLayer,fWire,fCell);
  
  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("-----------------------------------------------------------------------------------\n");
    printf("  Index  Time   Dt   Energy  PdgCode  PdgCode(M)  GenCode  SimID  McMom     \n");
    printf("-----------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;
  
  printf("%6i %10.3f %10.3f %10.5f %8i %8i %8i %8i %10.3f\n",
	 fIndex, 
	 fTime,
	 fDt,
	 fEnergy,
	 fPdgCode,
	 fMotherPdgCode,
	 fGeneratorCode,
	 fSimID,
	 fMcMomentum);
  
}
