//
#include <iostream>

#include "obj/TStnPid.hh"

ClassImp(TStnPid)

//-----------------------------------------------------------------------------
// Read an object of class TStnPid (version 1).
// V1 didn't store the cluster asymmetry data
//-----------------------------------------------------------------------------
// void TStnPid::ReadV1(TBuffer &R__b) {
// }

//-----------------------------------------------------------------------------
void TStnPid::Streamer(TBuffer& R__b) {
  int nwi, nwf;

  nwi = (int*  ) &fInt   - &fEleTrkNumber   + kNFreeInts;
  nwf = (float*) &fFloat - &fLogDedxProbEle + kNFreeFloats ;

  if (R__b.IsReading()) {
    //    Version_t R__v = R__b.ReadVersion(); 
    R__b.ReadVersion(); 
					// current version: V1
    R__b.ReadFastArray(&fEleTrkNumber  ,nwi);
    R__b.ReadFastArray(&fLogDedxProbEle,nwf);
  }
  else {
    R__b.WriteVersion(TStnPid::IsA());

    R__b.WriteFastArray(&fEleTrkNumber  ,nwi);
    R__b.WriteFastArray(&fLogDedxProbEle,nwf);
  }
}

//_____________________________________________________________________________
TStnPid::TStnPid(Int_t Number) {
  // 'Number' can be -1 ...
  
  fEleTrkNumber = Number;
}

//_____________________________________________________________________________
TStnPid::~TStnPid() {
}

//-----------------------------------------------------------------------------
void TStnPid::Clear(Option_t* opt) {
  Error("Print", "Not implemented yet");
}

//-----------------------------------------------------------------------------
void TStnPid::Print(Option_t* opt) const {
  Error("Print", "Not implemented yet");
}

// } // end namespace




