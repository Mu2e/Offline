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
// Read an object of class TStnCluster (version 1).
// V1 didn't store the cluster asymmetry data
//-----------------------------------------------------------------------------
void TStnCluster::ReadV1(TBuffer &R__b) {

  struct TStnClusterDataV1_t {
    int                       fNumber;             // index in the list of reconstructed clusters
    int                       fDiskID;	           // 
    int                       fNCrystals;          //
    int                       fNCr1     ;          // above 1 MeV
    int                       fTrackNumber;        // closest track in TStnTrackBlock
    int                       fIx1;	           // [row, column] or [x1,x2] for a disk
    int                       fIx2;
    int                       fInt[kNFreeIntsV1];
					           // floats
    float                     fX;
    float                     fY;
    float                     fZ;
    float                     fYMean;
    float                     fZMean;
    float                     fSigY;
    float                     fSigZ;
    float                     fSigR;
    float                     fEnergy;
    float                     fTime     ; 
    float                     fFrE1     ; // e1/etotal
    float                     fFrE2     ; // (e1+e2)/etotal
    float                     fSigE1    ;
    float                     fSigE2    ;
    float                     fFloat[kNFreeFloatsV1];
  };

  TStnClusterDataV1_t data;
  
  int                 nwi, nwf;
  
  nwi = (int*  ) data.fInt   - &data.fNumber + kNFreeIntsV1;
  nwf = (float*) data.fFloat - &data.fX      + kNFreeFloatsV1;
  
  R__b.ReadFastArray(&data.fNumber,nwi);
  R__b.ReadFastArray(&data.fX     ,nwf);
  
  fNumber        = data.fNumber      ;          // track index in the list of reconstructed clusters
  fDiskID        = data.fDiskID      ;	      // 
  fNCrystals     = data.fNCrystals   ;       //
  fNCr1          = data.fNCr1        ;       // above 1 MeV
  fTrackNumber   = data.fTrackNumber ;     // closest track in TStnTrackBlock
  fIx1           = data.fIx1         ;	      // [row, column] or [x1,x2] for a disk
  fIx2           = data.fIx2         ;
		  	       	          	// float part
  fX             = data.fX           ;
  fY             = data.fY           ;
  fZ             = data.fZ           ;
  fYMean         = data.fYMean       ;
  fZMean         = data.fZMean       ;
  fSigY          = data.fSigY        ; 
  fSigZ          = data.fSigZ        ;
  fSigR          = data.fSigR        ;
  fEnergy        = data.fEnergy      ;
  fTime          = data.fTime        ; 
  fFrE1          = data.fFrE1        ; // e1/etotal
  fFrE2          = data.fFrE2        ; // (e1+e2)/etotal
  fSigE1         = data.fSigE1       ;
  fSigE2         = data.fSigE2       ;
//-----------------------------------------------------------------------------
// initialize the V2 part. make sure teh numbers don't make sense
//-----------------------------------------------------------------------------
  fSigXX         = -1.;     // sums over crystals
  fSigXY         = -1.;
  fSigYY         = -1.;
  fNx            =  1.;     // cluster direction
  fNy            =  0.;
}

//-----------------------------------------------------------------------------
void TStnCluster::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = (int*  ) &fInt   - &fNumber + kNFreeInts;
  nwf = (float*) &fFloat - &fX      + kNFreeFloats ;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 

    if (R__v == 1) ReadV1(R__b);
    else {
					// current version: V2
      R__b.ReadFastArray(&fNumber,nwi);
      R__b.ReadFastArray(&fYMean ,nwf);
    }
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


//-----------------------------------------------------------------------------
void TStnCluster::Clear(Option_t* opt) {
  Error("Print", "Not implemented yet");
}

//-----------------------------------------------------------------------------
void TStnCluster::Print(Option_t* opt) const {
  Error("Print", "Not implemented yet");
}

// } // end namespace




