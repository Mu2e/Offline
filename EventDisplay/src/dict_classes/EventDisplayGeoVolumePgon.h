//
// Class which displays GeoVolumes with a hexagon (used e.g. by the Hexagon class). It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoPgon, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this vane.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_EventDisplayGeoVolumePgon_h
#define EventDisplay_src_dict_classes_EventDisplayGeoVolumePgon_h

#include <TGeoPgon.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "Offline/EventDisplay/src/EventDisplayFrame.h"
#include "Offline/EventDisplay/src/dict_classes/ComponentInfoContainer.h"
#include "Offline/EventDisplay/src/dict_classes/HistDraw.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class EventDisplayGeoVolumePgon : public TGeoVolume, public ComponentInfoContainer
{
  EventDisplayFrame *_mainframe;

  EventDisplayGeoVolumePgon();
  EventDisplayGeoVolumePgon(const EventDisplayGeoVolumePgon &);
  EventDisplayGeoVolumePgon& operator=(const EventDisplayGeoVolumePgon &);

  public:
#ifndef __CINT__
  EventDisplayGeoVolumePgon(double phi, double rmax, double halflength, EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfoContainer(info),_mainframe(mainframe)
  {
    //bare pointer needed since ROOT manages this object
    TGeoPgon *pgon=new TGeoPgon(nullptr, phi, 360.0, 6, 2);
    pgon->DefineSection(0,-halflength, 0, rmax);
    pgon->DefineSection(1,+halflength, 0, rmax);
    SetShape(pgon);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
  }
#endif

  virtual ~EventDisplayGeoVolumePgon()
  {
  }

  virtual const char* ClassName() const
  {
    TList  *l=IsA()->GetMenuList();
    l->Clear();

    TObject *obj = dynamic_cast<TObject*>(_mainframe);
    TClassMenuItem *m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),"Information","showInfo",obj,"TObject*",1); //TClassMenuItem accepts only bare pointers. m needs to be bare pointer because it is managed by ROOT.
    l->AddFirst(m);

    const std::vector<boost::shared_ptr<TObject> > &histVector = getComponentInfo()->getHistVector();
    unsigned int n = _mainframe->getHistDrawVector().size();
    for(unsigned int i=n; i<histVector.size(); i++) _mainframe->addHistDraw(); //add more HistDraws if there are more histograms in ComponentInfo

    const std::vector<boost::shared_ptr<HistDraw> > histDrawVector = _mainframe->getHistDrawVector();
    std::vector<boost::shared_ptr<TObject> >::const_iterator iter=histVector.begin();
    std::vector<boost::shared_ptr<HistDraw> >::const_iterator iter2=histDrawVector.begin();
    for(; iter!=histVector.end() && iter2!=histDrawVector.end(); iter++, iter2++)
    {
      m = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,IsA(),(*iter)->GetTitle(),"showHistogram",dynamic_cast<TObject*>(iter2->get()),"TObject*",1); //TClassMenuItem accepts only bare pointers. m needs to be bare pointer because it is managed by ROOT.
      l->Add(m);
    }

    IsA()->SetName(getComponentInfo()->getName()->c_str());
    return(IsA()->GetName());
  }

  ClassDef(EventDisplayGeoVolumePgon,0);
};

}
#endif /* EventDisplay_src_dict_classes_EventDisplayGeoVolumePgon_h */
