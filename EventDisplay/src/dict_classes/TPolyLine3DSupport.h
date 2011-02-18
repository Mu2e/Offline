//
// Class which displays lines of a tracker support structure (used, e.g. by the SupportTTracker class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this support structure. The context menu is overwritten with a menu item allowing the user to display information for this support structure.
//
// $Id: TPolyLine3DSupport.h,v 1.3 2011/02/18 04:10:55 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/18 04:10:55 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_SUPPORT_H
#define TPOLYLINE3D_SUPPORT_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DSupport : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DSupport();
  TPolyLine3DSupport(const TPolyLine3DSupport &);
  TPolyLine3DSupport& operator=(const TPolyLine3DSupport &);

  public:
#ifndef __CINT__
  TPolyLine3DSupport(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info)
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Support Structure and Envelope");
  }
#endif

  virtual ~TPolyLine3DSupport() {}

  ClassDef(TPolyLine3DSupport,0);
};

}
#endif
