//
// Class which builds the main frame for the event display, and provides functions to control the display, e.g. quit, moving to the next event, animations, storing the events into gif files (static and animated), detailed infos of tracks, hits, etc.
//
// $Id: EventDisplayFrame.h,v 1.6 2011/02/18 04:10:55 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/18 04:10:55 $
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
class TGRadioButton;
class TGLabel;
class TGTextEntry;
class TBox;
class TPolyLine;

namespace mu2e_eventdisplay
{
  class DataInterface;
  class EventDisplayPad;

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
    int            getEventToFind(bool &findEvent) const;
    void           showInfo(TObject*);
    void           fillZoomAngleFields();
    virtual void   CloseWindow(); //inherited from TGMainFrame 

    private:
    void updateHitLegend(bool draw);
    void updateTrackLegend(bool draw);
    void prepareAnimation();
    void drawSituation();
    void drawEverything();
    void combineAnimFiles();
    //the following functions are inherited from TGMainFrame
    virtual Bool_t HandleConfigureNotify(Event_t *event);
    virtual Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2);
    virtual Bool_t HandleTimer(TTimer *); //inherited function 
                                          //from TObject (gets called when 
                                          //timer times out - knows about 
                                          //this via TTimer::SetObject)

#ifndef __CINT__    //hide boost from ROOTCint
    boost::shared_ptr<DataInterface> _dataInterface;
#endif
    double              _timeCurrent, _timeStart, _timeStop;
    int                 _minHits, _eventToFind;
    bool                _isClosed, _findEvent;
    bool                _saveAnim;
    int                 _saveAnimCounter;
    std::string         _saveAnimFile;
    //bare pointers needed since ROOT manages these objects
    TRootEmbeddedCanvas *_mainCanvas, *_infoCanvas;
    EventDisplayPad     *_mainPad;
    TPad                *_infoPad;
    TText               *_clock;
    TTimer              *_timer;
    TGCheckButton       *_unhitButton, *_unhitCrystalsButton;
    TGCheckButton       *_supportStructuresButton, *_otherStructuresButton;
    TGCheckButton       *_outsideTracksButton, *_calorimeterViewButton, *_targetViewButton;
    TGCheckButton       *_hitColorButton, *_trackColorButton, *_backgroundButton;
    TGCheckButton       *_repeatAnimationButton;
    TGTextEntry         *_minHitField, *_eventToFindField;
    TGTextEntry         *_minXField, *_minYField, *_minZField, *_maxXField, *_maxYField, *_maxZField;
    TGTextEntry         *_phiField, *_thetaField, *_psiField;
    TGRadioButton       *_perspectiveButton, *_parallelButton;
    TGLabel             **_eventInfo;
    TText               *_legendText[30], *_legendParticleText[6];  
    TBox                *_legendBox[30];
    TPolyLine           *_legendParticleLine[6];  

    ClassDef(EventDisplayFrame,0);
  };
}

#endif

