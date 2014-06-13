//-----------------------------------------------------------------------------
// Mar 31 2000 P.Murat: minimal implementation of non-templated bitset class
//-----------------------------------------------------------------------------
#ifndef STNTUPLE_DATA_TBitset_hh
#define STNTUPLE_DATA_TBitset_hh

#include "TObject.h"

class TBitset: public TObject {
protected:
  Int_t     fNBits;
  Int_t     fNWords;		// 
  Int_t*    fBits;		// [fNWords]
public:
				// ****** constructors and destructor
  TBitset();
  TBitset(Int_t NBits);
  virtual ~TBitset();
				// ****** init methods

  Int_t         Init(Int_t NBits);

				// ****** accessors

  Int_t         GetNBits () const { return fNBits; }
  Int_t         GetNWords() const { return fNWords; }

  inline Int_t  GetBit(Int_t I) const ;

				// ****** modifiers

  inline void   SetBit (Int_t I);

  inline void   SetBit (Int_t I, int Val);  // Val = 0 or 1

				// ****** overloaded functions of TObject

  void  Clear(Option_t* Opt = "");

  ClassDef(TBitset,1)
};


//_____________________________________________________________________________
inline Int_t TBitset::GetBit(Int_t I) const {
  Int_t iw = I / 32;
  Int_t ib = I % 32;
  return (fBits[iw] >> ib) & 0x1 ;
}

//_____________________________________________________________________________
inline void TBitset::SetBit(Int_t I) {
  // this assumes that we are setting bits just once, certainly not safe,
  // I simply forgot C++ for XOR...

  Int_t iw = I / 32;
  Int_t ib = I % 32;
  fBits[iw] |= (0x1 << ib);
}

//_____________________________________________________________________________
inline void TBitset::SetBit(Int_t I, int Val) {
  // this assumes that we are setting bits just once, certainly not safe,
  // I simply forgot C++ for XOR...

  Val = (Val != 0);   // safety
  int iw = I / 32;
  int ib = I % 32;

  int mask = Val << ib;
  fBits[iw] = (fBits[iw]^(0x1 << ib)) | mask;
}


#endif

