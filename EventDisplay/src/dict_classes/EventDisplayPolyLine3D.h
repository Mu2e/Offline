//
// Class which displays 3D lines (used e.g. by Track, Cyliner class, etc.). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this track. The context menu is overwritten with a menu item allowing the user to display information for this track.
//
// $Id: EventDisplayPolyLine3D.h,v 1.1 2011/02/23 00:29:27 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/23 00:29:27 $
//
// Original author Ralf Ehrlich
//

#ifndef EVENTDISPLAY_POLYLINE3D_H
#define EVENTDISPLAY_POLYLINE3D_H

#include <TPolyLine3D.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class EventDisplayPolyLine3D : public TPolyLine3D, public ComponentInfo
{
  EventDisplayPolyLine3D();
  EventDisplayPolyLine3D(const EventDisplayPolyLine3D &);
  EventDisplayPolyLine3D& operator=(const EventDisplayPolyLine3D &);

  public:
#ifndef __CINT__
  EventDisplayPolyLine3D(const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TPolyLine3D(),ComponentInfo(info) 
  {
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
  }
#endif

  virtual ~EventDisplayPolyLine3D() {}

  virtual const char* ClassName() const
  {
    IsA()->SetName(getName()->c_str());
    return(IsA()->GetName());
  } 

  ClassDef(EventDisplayPolyLine3D,0);
};

}
#endif
