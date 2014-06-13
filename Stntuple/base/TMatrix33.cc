///////////////////////////////////////////////////////////////////////////////
// Apr 12 2001 P.Murat: fixed size matrix class to store track covariance 
// matrix. Constructor fakes dynamic memory allocation by TMatrix.
///////////////////////////////////////////////////////////////////////////////

#include "base/TMatrix33.hh"

ClassImp(TMatrix33)

//______________________________________________________________________________
void TMatrix33::Streamer(TBuffer &R__b)
{
   // Stream an object of class TMatrix33.

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TObject::Streamer(R__b);
    R__b.ReadFastArray(fElements,9);
  } 
  else {
    R__b.WriteVersion(TMatrix33::IsA());
    TObject::Streamer(R__b);
    R__b.WriteFastArray(fElements, 9);
  }
}

//_____________________________________________________________________________
TMatrix33::TMatrix33() {
  // allocate new matrix. Arguments are number of rows, columns, row
  // lowerbound (0 default) and column lowerbound (0 default). Fake memory
  // allocation by TMatrix class and as such use static memory alocation

  Allocate(3,3,0,0,1);
}


//_____________________________________________________________________________
TMatrix33::~TMatrix33() {
  // destructor: does nothing ! (has to overload destructor of TMatrix)
  Invalidate();
}


