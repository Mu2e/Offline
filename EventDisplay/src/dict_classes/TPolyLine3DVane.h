//
// Class which displays a vane of the calometer (used by Vane class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this vane. The context menu is overwritten with a menu item allowing the user to display information for this vane.
//
// $Id: TPolyLine3DVane.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_VANE_H
#define TPOLYLINE3D_VANE_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DVane : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DVane();
  TPolyLine3DVane(const TPolyLine3DVane &);
  TPolyLine3DVane& operator=(const TPolyLine3DVane &);

  public:
#ifndef __CINT__
  TPolyLine3DVane(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info)
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Vane");
  }
#endif

  virtual ~TPolyLine3DVane() {}

  ClassDef(TPolyLine3DVane,0);
};

}
#endif
