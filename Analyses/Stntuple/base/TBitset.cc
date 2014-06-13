//-----------------------------------------------------------------------------
// Mar 31 2000 P.Murat: minimal implementation of non-templated bitset class
// - Jul 16 2001: for nbits=32,64 etc one extra word will be allocated,
//                as the plan is to replace it with the ROOT 3.01 class, 
//                leave it as it is for the time being
//-----------------------------------------------------------------------------
#include "base/TBitset.hh"

ClassImp(TBitset)

//______________________________________________________________________________
void TBitset::Streamer(TBuffer &R__b)
{
  // Stream an object of class TBitset, don't write out fNWords, because it is
  // a derivative from fNBits. Also dont' write out TObject part.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      R__b >> fNBits;
      if (fNBits == 0) {
	if (fBits) {
	  delete [] fBits; 
	  fBits = 0;
	}
	fNWords = 0;
      } else {
	int nw  = fNBits/32+1;
	if (nw != fNWords) {
	  delete [] fBits; 
	  fNWords = nw;
	  fBits   = new Int_t[fNWords]; 
	}
      }
      R__b.ReadFastArray(fBits,fNWords); 
   } else {
      R__b.WriteVersion(TBitset::IsA());
      R__b << fNBits;
      R__b.WriteFastArray(fBits,fNWords); 
   }
}

//_____________________________________________________________________________
TBitset::TBitset() {
  fNBits  = 0;
  fNWords = 0;
  fBits   = NULL;
}

//_____________________________________________________________________________
TBitset::TBitset(Int_t NBits): 
  fNBits(NBits),
  fNWords(NBits/32+1)
{
  // by default set bitset contents to zero
  fBits = new Int_t[fNWords];
  memset(fBits,0,fNWords*4);
}

//_____________________________________________________________________________
Int_t TBitset::Init(Int_t NBits)
{
  // by default set bitset contents to zero

  if (fNBits != NBits) {
    
    if (fBits != 0) delete [] fBits;
    fNBits = NBits;
    fNWords = fNBits/32+1;
    fBits   = new Int_t [fNWords];
  }
  memset(fBits,0,fNWords*4);
  return 0;
}

//_____________________________________________________________________________
TBitset::~TBitset() {
  if (fBits) delete [] fBits;
}

//_____________________________________________________________________________
void TBitset::Clear(Option_t* Opt) {
  memset(fBits,0,fNWords*4);
}

