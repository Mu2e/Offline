//
// Class which displays a cylinder (used, e.g. by the Cyliner class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoTube, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this support structure.
//
// $Id: EventDisplayGeoVolumeTube.h,v 1.4 2011/09/04 04:43:34 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/09/04 04:43:34 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_EventDisplayGeoVolumeTube_h
#define EventDisplay_src_dict_classes_EventDisplayGeoVolumeTube_h

#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class EventDisplayGeoVolumeTube : public TGeoVolume, public ComponentInfo
{
  public:
  EventDisplayGeoVolumeTube();
  EventDisplayGeoVolumeTube(const EventDisplayGeoVolumeTube &);
  EventDisplayGeoVolumeTube& operator=(const EventDisplayGeoVolumeTube &);

  public:
#ifndef __CINT__
  EventDisplayGeoVolumeTube(double innerRadius, double outerRadius, double halflength, const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfo(info)
  {
    //bare pointer needed since ROOT manages this object
    TGeoTube *tube=new TGeoTube(NULL, innerRadius, outerRadius, halflength);
    SetShape(tube);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
  }
#endif

  virtual ~EventDisplayGeoVolumeTube()
  {
  }

  virtual const char* ClassName() const
  {
    IsA()->SetName(getName()->c_str());
    return(IsA()->GetName());
  }

  ClassDef(EventDisplayGeoVolumeTube,0);
};

}
#endif /* EventDisplay_src_dict_classes_EventDisplayGeoVolumeTube_h */
