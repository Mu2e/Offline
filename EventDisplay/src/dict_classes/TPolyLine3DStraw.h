//
// Class which displays straws (used by Straw class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this straw/hit. The context menu is overwritten with a menu item allowing the user to display information for this straw/hit.
//
// $Id: TPolyLine3DStraw.h,v 1.1 2011/01/26 18:11:44 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/01/26 18:11:44 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_STRAW_H
#define TPOLYLINE3D_STRAW_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DStraw : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DStraw();
  TPolyLine3DStraw(const TPolyLine3DStraw &);
  TPolyLine3DStraw& operator=(const TPolyLine3DStraw &);

  public:
#ifndef __CINT__
  TPolyLine3DStraw(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info)
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Straw");
  }
#endif

  virtual ~TPolyLine3DStraw() {}

  ClassDef(TPolyLine3DStraw,0);
};

}
#endif
