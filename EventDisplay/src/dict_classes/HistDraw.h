//
// Class which displays histograms and graphs. The reason why this class is necessary, is because there will be a variable number of histograms, which all need the same function to display them (e.g. the ShowHistogram function). However, TClassMenuItem can be associated only with a function with one parameter. The parameter needs to be clicked object (e.g. the polyline, etc.) which contains the ContentInfo, which holds ALL histograms. Since one cannot pass a parameter to this function to tell it which of the histograms it needs to drawn, one has to use a more complex approach: The functions (showHistogram) which are used for TClassMenuItems belong to special object of the class HistDraw. Each of these objects has a unique number (_index) which tell this function which of the histogram it needs to draw. These objects are given to the TClassMenuItem. The number of such objects depends on the number of histograms. The objects are stored in a vector which is owned by the EventDisplayFrame object.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_HistDraw_h
#define EventDisplay_src_dict_classes_HistDraw_h

#include <TPad.h>
#include <TRootEmbeddedCanvas.h>
#include <TGCanvas.h>
#include "EventDisplay/src/dict_classes/ComponentInfoContainer.h"

namespace mu2e_eventdisplay
{

class HistDraw : public TObject  //needs to inherit from TObject to work with TClassMenuItem
{
  int                 _index;
    //bare pointers needed since ROOT manages these objects
  TPad                *_infoPad, *_mainPad;
  TRootEmbeddedCanvas *_infoEmbeddedCanvas;
  TGCanvas            *_infoCanvas;
  
  HistDraw();
  HistDraw(const HistDraw &);
  HistDraw& operator=(const HistDraw &);

  public:
  HistDraw(int index, TPad *infoPad, TPad *mainPad, TRootEmbeddedCanvas *infoEmbeddedCanvas, TGCanvas *infoCanvas) : 
           _index(index), _infoPad(infoPad), _mainPad(mainPad), 
           _infoEmbeddedCanvas(infoEmbeddedCanvas), _infoCanvas(infoCanvas) {}
  void showHistogram(TObject *o)
  {
    ComponentInfoContainer *container=dynamic_cast<ComponentInfoContainer*>(o);
    if(!container) return;

    _infoPad->cd();
    _infoPad->Clear();
    unsigned int width=_infoCanvas->GetWidth()-20;
    unsigned int height=_infoCanvas->GetHeight()-20;
    _infoEmbeddedCanvas->SetWidth(width);
    _infoEmbeddedCanvas->SetHeight(height);
    _infoCanvas->Layout();
    _infoPad->cd();
    _infoPad->Range(0, 0, 1, 1);
    _infoPad->Modified();
    _infoPad->Update();
    container->getComponentInfo()->showHist(_index);
    _infoPad->Modified();
    _infoPad->Update();
    _mainPad->cd();
  }

  ClassDef(HistDraw,0);
};

}

#endif
