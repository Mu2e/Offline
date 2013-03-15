//
// Class which displays a Cone. It is inherited from ROOT's TGeoVolume and the ComponentInfo class which stores specific information for this support structure. The class' constructure creates a TGeoCone, which is put into the TGeoVolume. The context menu is overwritten with a menu item allowing the user to display information for this support structure.
//
// $Id: EventDisplayGeoVolumeCone.h,v 1.2 2013/03/15 16:20:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 16:20:00 $
//
// Original author MyeongJae Lee, based on Ralf Ehrlich's EventDisplayGeoVolumeTube.h codes
//

#ifndef EventDisplay_src_dict_classes_EventDisplayGeoVolumeCone_h
#define EventDisplay_src_dict_classes_EventDisplayGeoVolumeCone_h

#include <TGeoCone.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include "EventDisplay/src/EventDisplayFrame.h"
#include "EventDisplay/src/dict_classes/ComponentInfoContainer.h"
#include "EventDisplay/src/dict_classes/HistDraw.h"
#include <TClass.h>
#include <TList.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class EventDisplayGeoVolumeCone : public TGeoVolume, public ComponentInfoContainer
{
  EventDisplayFrame *_mainframe;

  EventDisplayGeoVolumeCone();
  EventDisplayGeoVolumeCone(const EventDisplayGeoVolumeCone &);
  EventDisplayGeoVolumeCone& operator=(const EventDisplayGeoVolumeCone &);

  public:
#ifndef __CINT__
  EventDisplayGeoVolumeCone(double innerRadius1, double outerRadius1, double innerRadius2, double outerRadius2, double halflength, EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info):TGeoVolume(),ComponentInfoContainer(info),_mainframe(mainframe)
  {
    //bare pointer needed since ROOT manages this object
    TGeoCone *cone=new TGeoCone(nullptr, halflength, innerRadius1, outerRadius1, innerRadius2, outerRadius2);
    SetShape(cone);
    SetNumber(GetGeoManager()->AddVolume(this)); //this is what happens in the base TGeoVolume constructor
  }
#endif

  virtual ~EventDisplayGeoVolumeCone()
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

  ClassDef(EventDisplayGeoVolumeCone,0);
};

}
#endif /* EventDisplay_src_dict_classes_EventDisplayGeoVolumeCone_h */
