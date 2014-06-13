// $Id: TNamedHandle.hh,v 1.1 2014/06/13 06:13:04 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:13:04 $
//
// Contact person Pavel Murat
//
#ifndef murat_inc_TNamedHandle_hh
#define murat_inc_TNamedHandle_hh


#include "TNamed.h"

class TNamedHandle : public TNamed {
  public:
  void*    fObject;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TNamedHandle();
  TNamedHandle(const char* Name, void* Object);
  ~TNamedHandle();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void* Object() { return fObject; }

  ClassDef(TNamedHandle,0)

};

#endif
