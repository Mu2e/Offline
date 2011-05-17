//
// Class which displays 3D lines (used e.g. by Track, Cyliner class, etc.). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this track. The context menu is overwritten with a menu item allowing the user to display information for this track.
//
// $Id: EventDisplayPolyLine3D.h,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:35 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h
#define EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h

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
#endif /* EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h */
