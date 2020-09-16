//
// Class which displays 3D lines (used e.g. by Track, Cyliner class, etc.). It is inherited from ROOT's TPolyLine3D and the ComponentInfo class which stores specific information for this track. The context menu is overwritten with a menu item allowing the user to display information for this track.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h
#define EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h

#include <TPolyLine3D.h>
#include "EventDisplay/src/EventDisplayFrame.h"
#include "EventDisplay/src/dict_classes/ComponentInfoContainer.h"
#include "EventDisplay/src/dict_classes/HistDraw.h"
#include <TClass.h>
#include <TClassMenuItem.h>

namespace mu2e_eventdisplay
{

class EventDisplayPolyLine3D : public TPolyLine3D, public ComponentInfoContainer
{
  EventDisplayFrame *_mainframe;

  EventDisplayPolyLine3D();
  EventDisplayPolyLine3D(const EventDisplayPolyLine3D &);
  EventDisplayPolyLine3D& operator=(const EventDisplayPolyLine3D &);

  public:
#ifndef __CINT__
  EventDisplayPolyLine3D(EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info):
                          TPolyLine3D(),ComponentInfoContainer(info),_mainframe(mainframe)
  {}
#endif

  virtual ~EventDisplayPolyLine3D() {}

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

  ClassDef(EventDisplayPolyLine3D,0);
};

}
#endif /* EventDisplay_src_dict_classes_EventDisplayPolyLine3D_h */
