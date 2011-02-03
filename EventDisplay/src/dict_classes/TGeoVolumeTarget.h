//
// Class which displays the cylinder of the target (used by the Target class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoTube, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this target.
//
// $Id: TGeoVolumeTarget.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef TGEOTUBE_TARGET_H
#define TGEOTUBE_TARGET_H

#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TGeoVolumeTarget : public TGeoVolume, public ComponentInfo
{
  TGeoVolumeTarget();
  TGeoVolumeTarget(const TGeoVolumeTarget &);
  TGeoVolumeTarget& operator=(const TGeoVolumeTarget &);

  public:
#ifndef __CINT__
  TGeoVolumeTarget(double innerRadius, double outerRadius, double halflength, const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfo(info)
  {
    //bare pointer needed since ROOT manages this object
    TGeoTube *tube=new TGeoTube(NULL, innerRadius, outerRadius, halflength);
    SetShape(tube);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Target");
  }
#endif

  virtual ~TGeoVolumeTarget() 
  {
  }

  ClassDef(TGeoVolumeTarget,0);
};

}
#endif
