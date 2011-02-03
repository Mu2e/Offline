//
// Class which displays a crystal of the calometer (used by Crystal class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this crystal. The context menu is overwritten with a menu item allowing the user to display information for this crystal.
//
// $Id: TPolyLine3DCrystal.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_CRYSTAL_H
#define TPOLYLINE3D_CRYSTAL_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DCrystal : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DCrystal();
  TPolyLine3DCrystal(const TPolyLine3DCrystal &);
  TPolyLine3DCrystal& operator=(const TPolyLine3DCrystal &);

  public:
#ifndef __CINT__
  TPolyLine3DCrystal(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info)
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Crystal");
  }
#endif

  virtual ~TPolyLine3DCrystal() {}

  ClassDef(TPolyLine3DCrystal,0);
};

}
#endif
