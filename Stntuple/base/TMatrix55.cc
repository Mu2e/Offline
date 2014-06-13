///////////////////////////////////////////////////////////////////////////////
// Apr 12 2001 P.Murat: fixed size matrix class to store track covariance 
// matrix. Constructor fakes dynamic memory allocation by TMatrix.
///////////////////////////////////////////////////////////////////////////////

#include "base/TMatrix55.hh"

ClassImp(TMatrix55)

//______________________________________________________________________________
void TMatrix55::Streamer(TBuffer &R__b)
{
   // Stream an object of class TMatrix55.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TObject::Streamer(R__b);
    R__b.ReadFastArray(fElements,25);
  } 
  else {
    R__b.WriteVersion(TMatrix55::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(fElements,25);
  }
}

//_____________________________________________________________________________
TMatrix55::TMatrix55() {
  // allocate new matrix. Arguments are number of rows, columns, row
  // lowerbound (0 default) and column lowerbound (0 default). Fake memory
  // allocation by TMatrix class and as such use static memory alocation

  Allocate(5,5,0,0,1);
}


//_____________________________________________________________________________
TMatrix55::~TMatrix55() {
  // destructor: does nothing ! (has to overload destructor of TMatrix)
  Invalidate();
}


