//
// Class which displays the cylinder of a tracker support structure (used, e.g. by the SupportTTracker class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoTube, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this support structure.
//
// $Id: TGeoVolumeSupport.h,v 1.3 2011/02/18 04:10:55 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/18 04:10:55 $
//
// Original author Ralf Ehrlich
//

#ifndef TGEOTUBE_SUPPORT_H
#define TGEOTUBE_SUPPORT_H

#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "ComponentInfo.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class TGeoVolumeSupport : public TGeoVolume, public ComponentInfo
{
  TGeoVolumeSupport();
  TGeoVolumeSupport(const TGeoVolumeSupport &);
  TGeoVolumeSupport& operator=(const TGeoVolumeSupport &);

  public:
#ifndef __CINT__
  TGeoVolumeSupport(double innerRadius, double outerRadius, double halflength, const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfo(info)
  {
    //bare pointer needed since ROOT manages this object
    TGeoTube *tube=new TGeoTube(NULL, innerRadius, outerRadius, halflength);
    SetShape(tube);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
    TList  *l=IsA()->GetMenuList();
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",const_cast<TObject*>(mainframe),"TObject*",1); //bare pointer needed since ROOT manages this object
    l->Clear();
    l->AddFirst(m);
    IsA()->SetName("Support Structure and Envelope");
  }
#endif

  virtual ~TGeoVolumeSupport() 
  {
  }

  ClassDef(TGeoVolumeSupport,0);
};

}
#endif
