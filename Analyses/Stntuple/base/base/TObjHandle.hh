// $Id: TObjHandle.hh,v 1.1 2014/06/13 05:24:57 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 05:24:57 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_base_TObjHandle_hh
#define Stntuple_base_TObjHandle_hh


#include "TObject.h"

class TObjHandle : public TObject {
  public:
  void*    fObject;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TObjHandle();
  TObjHandle(void* Object);
  ~TObjHandle();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void* Object() { return fObject; }

  ClassDef(TObjHandle,0)

};

#endif
