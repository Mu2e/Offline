#ifndef TMatrix55_hh
#define TMatrix55_hh
//-----------------------------------------------------------------------------
//  Apr 18 2001 P.Murat: need fixed size matrix class to store track 
//  cov. matrix which wouldn't be allocating memory dynamically
//-----------------------------------------------------------------------------

#include "TMatrixFSym.h"
class TMatrix55 : public TMatrixFSym {
public:
				// constructor and destructor, all the rest
				// methods are coming from TMatrix...
  TMatrix55();
  virtual ~TMatrix55();

  ClassDef(TMatrix55,3)
};

#endif
