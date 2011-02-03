//
// Class which displays a vane of the calometer (used by the Vane class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoBox, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this vane.
//
// $Id: TGeoVolumeVane.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef TGEOTUBE_VANE_H
#define TGEOTUBE_VANE_H

#include <TGeoBBox.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TGeoVolumeVane : public TGeoVolume, public ComponentInfo
{
  TGeoVolumeVane();
  TGeoVolumeVane(const TGeoVolumeVane &);
  TGeoVolumeVane& operator=(const TGeoVolumeVane &);

  public:
#ifndef __CINT__
  TGeoVolumeVane(double dx, double dy, double dz, const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfo(info)
  {
    //bare pointer needed since ROOT manages this object
    TGeoBBox *box=new TGeoBBox(NULL, dx, dy, dz);
    SetShape(box);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Vane");
  }
#endif

  virtual ~TGeoVolumeVane() 
  {
  }

  ClassDef(TGeoVolumeVane,0);
};

}
#endif
