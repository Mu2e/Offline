//
// Class which displays a steelshield (used by the SteelShield class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoBox, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this vane.
//
// $Id: TGeoVolumeSteelShield.h,v 1.1 2011/02/18 04:10:55 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/18 04:10:55 $
//
// Original author Ralf Ehrlich
//

#ifndef TGEOTUBE_STEELSHIELD_H
#define TGEOTUBE_STEELSHIELD_H

#include <TGeoBBox.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TGeoVolumeSteelShield : public TGeoVolume, public ComponentInfo
{
  TGeoVolumeSteelShield();
  TGeoVolumeSteelShield(const TGeoVolumeSteelShield &);
  TGeoVolumeSteelShield& operator=(const TGeoVolumeSteelShield &);

  public:
#ifndef __CINT__
  TGeoVolumeSteelShield(double dx, double dy, double dz, const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfo(info)
  {
    //bare pointer needed since ROOT manages this object
    TGeoBBox *box=new TGeoBBox(NULL, dx, dy, dz);
    SetShape(box);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("CR Steel Shield");
  }
#endif

  virtual ~TGeoVolumeSteelShield() 
  {
  }

  ClassDef(TGeoVolumeSteelShield,0);
};

}
#endif
