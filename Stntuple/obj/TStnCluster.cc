//
#include <iostream>

#include "TMatrix.h"
#include "obj/TStnCluster.hh"

// #include "RecoDataProducts/inc/CaloCluster.hh"

namespace {
  //  const double BF = 1.4116 ;  // CDF case
  const double BF = 1.0 ;         // MU2E case
}

//namespace murat {

ClassImp(TStnCluster)

//-----------------------------------------------------------------------------
void TStnCluster::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = ((int*  ) &fYMean )      - &fNumber;
  nwf = ((float*) &fCaloCluster) - &fYMean ;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }

    R__b.ReadFastArray(&fNumber,nwi);
    R__b.ReadFastArray(&fYMean ,nwf);
  }
  else {
    R__b.WriteVersion(TStnCluster::IsA());

    R__b.WriteFastArray(&fNumber,nwi);
    R__b.WriteFastArray(&fYMean ,nwf);
  }
}

//_____________________________________________________________________________
TStnCluster::TStnCluster(Int_t Number) {
  // 'Number' can be -1 ...
  
  fNumber       = Number;
  fCaloCluster  = 0;
  fClosestTrack = 0;
}


//_____________________________________________________________________________
TStnCluster::~TStnCluster() {
}


//_____________________________________________________________________________
void TStnCluster::Print(Option_t* opt) const {
  Error("Print", "Not implemented yet");
}

// } // end namespace




