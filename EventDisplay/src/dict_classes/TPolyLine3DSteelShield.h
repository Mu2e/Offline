//
// Class which displays a steel shield (used by SteelShield class). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this vane. The context menu is overwritten with a menu item allowing the user to display information for this vane.
//
// $Id: TPolyLine3DSteelShield.h,v 1.1 2011/02/18 04:10:55 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/18 04:10:55 $
//
// Original author Ralf Ehrlich
//

#ifndef TPOLYLINE3D_STEELSHIELD_H
#define TPOLYLINE3D_STEELSHIELD_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TPolyLine3DSteelShield : public TPolyLine3D, public ComponentInfo
{
  TPolyLine3DSteelShield();
  TPolyLine3DSteelShield(const TPolyLine3DSteelShield &);
  TPolyLine3DSteelShield& operator=(const TPolyLine3DSteelShield &);

  public:
#ifndef __CINT__
  TPolyLine3DSteelShield(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info)
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("CR Steel Shield");
  }
#endif

  virtual ~TPolyLine3DSteelShield() {}

  ClassDef(TPolyLine3DSteelShield,0);
};

}
#endif
