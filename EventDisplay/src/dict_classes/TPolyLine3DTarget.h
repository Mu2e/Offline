//
// Class which displays the target (used by Target class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this target. The context menu is overwritten with a menu item allowing the user to display information for this target.
//
// $Id: TPolyLine3DTarget.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_TARGET_H
#define TPOLYLINE3D_TARGET_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DTarget : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DTarget();
  TPolyLine3DTarget(const TPolyLine3DTarget &);
  TPolyLine3DTarget& operator=(const TPolyLine3DTarget &);

  public:
#ifndef __CINT__
  TPolyLine3DTarget(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info) 
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Target");
  }
#endif

  virtual ~TPolyLine3DTarget() {}

  ClassDef(TPolyLine3DTarget,0);
};

}
#endif
