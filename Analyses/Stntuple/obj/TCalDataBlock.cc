///////////////////////////////////////////////////////////////////////////////
//  Dec 07 2001 P.Murat: start putting in some comments
//  ---------------------------------------------------
// TCalDataBlock: ROOT-parseable description of TCalData to be stored in 
//                STNTUPLE
///////////////////////////////////////////////////////////////////////////////
#include "TVector2.h"

#include "Stntuple/obj/TCalHitData.hh"
#include "Stntuple/obj/TCalDataBlock.hh"

ClassImp(TCalDataBlock)

//_____________________________________________________________________________
void TCalDataBlock::ReadV1(TBuffer &R__b) {

  struct TCalDataBlockV1_t {
    int            fNHits;		// number of hit crystals
    int            fNDisks;             // 
    int            fNCrystals  [4];	// 
    float          fRMin       [4];	// as a temporary measure, store 
    float          fRMax       [4];
    float          fZ0         [4];     // 
    float          fCrystalSize;
    float          fMinFraction;        // min fr of the included crystal area

    TClonesArray*  fListOfCalHitData;	// list of crystal hit data 
  };

  TCalDataBlockV1_t data; 

  int nwi = ((int*  ) data.fRMin       ) - &data.fNHits;
  int nwf = ((float*) &data.fListOfCalHitData) - data.fRMin;

  R__b.ReadFastArray(&fNHits,nwi);
  R__b.ReadFastArray(fRMin  ,nwf);

  if (fNHits > 0) {
    fListOfCalHitData->Streamer(R__b);
  }
				// initialize V2 variables 
  fWrapperThickness = 0.065;    // 65 microns
  fShellThickness   = 0;
}

//_____________________________________________________________________________
void TCalDataBlock::ReadV2(TBuffer &R__b) {

  struct TCalDataBlockV2_t {
    int            fNHits;		// number of hit crystals
    int            fNDisks;             // 
    int            fNCrystals  [4];	// 
    float          fRMin       [4];	// as a temporary measure, store 
    float          fRMax       [4];
    float          fZ0         [4];     // 
    float          fCrystalSize;
    float          fMinFraction;        // min fr of the included crystal area
    float          fWrapperThickness;   //
    float          fShellThickness;     //

    TClonesArray*  fListOfCalHitData;	// list of crystal hit data 
  };

  TCalDataBlockV2_t data; 

  int nwi = ((int*  ) data.fRMin             ) - &data.fNHits;
  int nwf = ((float*) &data.fListOfCalHitData) - data.fRMin;

  R__b.ReadFastArray(&fNHits,nwi);
  R__b.ReadFastArray(fRMin  ,nwf);

  if (fNHits > 0) {
    fListOfCalHitData->Streamer(R__b);
  }

}

//______________________________________________________________________________
void TCalDataBlock::Streamer(TBuffer &R__b)
{
  // Stream an object of class TCalDataBlock.
  // V2 adds out of time energies
  // V3 adds miniplug variables (code by Angela Wyatt)

  int nwi = ((int*  ) fRMin             ) - &fNHits;
  int nwf = ((float*) &fListOfCalHitData) - fRMin;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    if      (R__v == 1) ReadV1(R__b);
    else if (R__v == 2) ReadV2(R__b);
    else {
//-----------------------------------------------------------------------------
// read version > 1 ???
//-----------------------------------------------------------------------------
      printf(">>> ERROR: TCalDataBlock::Streamer read wrong version %i\n",R__v);
    } 
  }
  else {
    R__b.WriteVersion(TCalDataBlock::IsA());

    R__b.WriteFastArray(&fNHits,nwi);
    R__b.WriteFastArray(fRMin  ,nwf);

    if (fNHits > 0) {
      fListOfCalHitData->Streamer(R__b);
    }
  }
}

//_____________________________________________________________________________
TCalDataBlock::TCalDataBlock() {
				// commented out pieces are available starting 
				// from version 2.24/05

  fListOfCalHitData    = new TClonesArray("TCalHitData",100);
  fListOfCalHitData->BypassStreamer(kFALSE);
  //  fAdcThreshold = 0;
  Clear();
}

//_____________________________________________________________________________
TCalDataBlock::~TCalDataBlock() {
  fListOfCalHitData->Delete();
  delete fListOfCalHitData;
}


//_____________________________________________________________________________
void TCalDataBlock::Clear(Option_t* opt) {
  fListOfCalHitData->Clear();
  fNHits      = 0;
}


//_____________________________________________________________________________
void TCalDataBlock::Print(Option_t* opt) const {
  // print all the towers in the list
  if (fNHits > 0) {
    fListOfCalHitData->At(0)->Print("banner");
    for (int i=0; i<fNHits; i++) {
      fListOfCalHitData->At(i)->Print();
    }
  }
}


//-----------------------------------------------------------------------------
int TCalDataBlock::DiskNumber(int CrystalID) {
  printf(">>> ERROR TCalDataBlock::DiskNumber NOT IMPLEMENTED YET\n");
  return -1;
}

//-----------------------------------------------------------------------------
int TCalDataBlock::CrystalRadius(int CrystalID) {
  printf(">>> ERROR TCalDataBlock::CrystalRadius NOT IMPLEMENTED YET\n");
  return -1;
}
