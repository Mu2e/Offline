//
// Class which builds the main frame for the event display, and provides functions to control the display, e.g. quit, moving to the next event, animations, storing the events into gif files (static and animated), detailed infos of tracks, hits, etc.
//
// $Id: EventDisplayFrame.h,v 1.2 2011/01/29 02:14:20 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/01/29 02:14:20 $
//
// Original author Ralf Ehrlich
//

#ifndef MU2EMAINFRAME_H
#define MU2EMAINFRAME_H

#include <iostream>
#include <TGFrame.h>
#ifndef __CINT__
#include "FWCore/Framework/interface/Event.h"
#include "boost/shared_ptr.hpp"
#endif

class TRootEmbeddedCanvas;
class TPad;
class TText;
class TTimer;
class TGCheckButton;
class TGLabel;
class TGTextEntry;
class TBox;
class TPolyLine;

namespace mu2e_eventdisplay
{
  class DataInterface;

  class EventDisplayFrame : public TGMainFrame
  {
    EventDisplayFrame();
    EventDisplayFrame(const EventDisplayFrame &);
    EventDisplayFrame& operator=(const EventDisplayFrame &);

    public:
    EventDisplayFrame(const TGWindow* p, UInt_t w, UInt_t h);
    virtual        ~EventDisplayFrame();
    void           fillGeometry();
#ifndef __CINT__   //hide edm::Event from ROOTCint
    void           fillEvent(const edm::Event& event);
#endif
    bool           isClosed() const;
    int            getMinimumHits() const;
    void           prepareAnimation();
    void           showInfo(TObject*);
    //the following functions are inherited from TGMainFrame
    Bool_t         HandleConfigureNotify(Event_t *event);
    virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);
    virtual void   CloseWindow();  
    Bool_t         HandleTimer(TTimer *); //inherited function 
                                          //from TObject (gets called when 
                                          //timer times out - knows about 
                                          //this via TTimer::SetObject)

    private:
    void updateHitLegend(bool draw);
    void updateTrackLegend(bool draw);
    void drawSituation();
    void drawEverything();
    void combineAnimFiles();

#ifndef __CINT__    //hide boost from ROOTCint
    boost::shared_ptr<DataInterface> _dataInterface;
#endif
    double              _timeCurrent, _timeStart, _timeStop;
    int                 _minHits;
    bool                _isClosed;
    bool                _saveAnim;
    int                 _saveAnimCounter;
    std::string         _saveAnimFile;
    //bare pointers needed since ROOT manages these objects
    TRootEmbeddedCanvas *_mainCanvas, *_infoCanvas;
    TPad                *_mainPad, *_infoPad;
    TText               *_clock;
    TTimer              *_timer;
    TGCheckButton       *_unhitButton, *_supportStructureButton, *_outsideTracksButton;
    TGCheckButton       *_hitColorButton, *_trackColorButton, *_backgroundButton;
    TGTextEntry         *_minHitField;
    TGLabel             **_eventInfo;
    TText               *_legendText[30], *_legendParticleText[6];  
    TBox                *_legendBox[30];
    TPolyLine           *_legendParticleLine[6];  

    ClassDef(EventDisplayFrame,0);
  };
}

#endif

