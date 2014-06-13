#ifndef TMatrix33_hh
#define TMatrix33_hh
//-----------------------------------------------------------------------------
//  Apr 18 2001 P.Murat: need fixed size matrix class to store 3x3
//  cov. matrices which wouldn't allocate memory dynamically
//-----------------------------------------------------------------------------

#include "TMatrixFSym.h"
class TMatrix33 : public TMatrixFSym {
public:
				// constructor and destructor, all the rest
				// methods are coming from TMatrix...
  TMatrix33();
  virtual ~TMatrix33();

  ClassDef(TMatrix33,3)
};

#endif
