// $Id: TObjHandle.hh,v 1.1 2014/06/13 06:13:04 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:13:04 $
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
