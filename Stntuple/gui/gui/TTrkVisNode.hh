#ifndef TTrkVisNode_hh
#define TTrkVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TArc.h"
//------------------------------------------------------------------------------
// this clause is necessary
//-----------------------------------------------------------------------------
#ifndef __CINT__
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#else
namespace mu2e {
  class StrawHitCollection;
  class StrawHitPositionCollection;
  class StrawHitFlagCollection;
  class PtrStepPointMCVectorCollection;
  class StrawDigiMCCollection;
  class CalTimePeakCollection;
  class CalTimePeak;
  class KalRepPtrCollection;
  class TTracker;
};
#endif

#include "Stntuple/base/TVisNode.hh"

class TStnTrackBlock;
class TEvdStrawTracker;

class TTrkVisNode: public TVisNode {
public:
  enum {
    kPickHits     = 0,
    kPickTracks   = 1,
    kPickClusters = 2
  };
  
protected:

  const mu2e::StrawHitCollection**             fStrawHitColl;
  const mu2e::StrawHitPositionCollection**     fStrawHitPosColl;  //
  const mu2e::StrawHitFlagCollection**         fStrawHitFlagColl; //
  const mu2e::CalTimePeakCollection**          fCalTimePeakColl;  //
  const mu2e::PtrStepPointMCVectorCollection** fMcPtrColl; 
  const mu2e::StrawDigiMCCollection**          fStrawDigiMCColl; 
  const mu2e::KalRepPtrCollection**            fKalRepPtrColl;

  TStnTrackBlock*   fTrackBlock;
  Color_t           fTrackColor;

  TEvdStrawTracker* fTracker;

  TArc*                     fArc;
  const mu2e::CalTimePeak*  fTimePeak;

  Int_t         fDisplayBackgroundHits;
  Int_t         fTimeWindow;
  Int_t         fPickMode;
  Int_t         fUseStereoHits;
  double        fEventTime;

  TObjArray*    fListOfStrawHits;
  TObjArray*    fListOfTracks;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TTrkVisNode();
  TTrkVisNode(const char* Name, const mu2e::TTracker* Tracker, TStnTrackBlock* fTrackBlock);
  virtual ~TTrkVisNode();

//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObjArray* GetListOfTracks() { return fListOfTracks; }
  Color_t    GetTrackColor  () { return fTrackColor;   }

  const mu2e::StrawHitCollection* GetStrawHitColl() { 
    return *fStrawHitColl; 
  }

  const mu2e::StrawHitPositionCollection* GetStrawHitPosColl() { 
    return *fStrawHitPosColl;
  }

  const mu2e::StrawHitFlagCollection* GetStrawHitFlagColl() { 
    return *fStrawHitFlagColl;
  }

  const mu2e::PtrStepPointMCVectorCollection* GetMcPtrColl() { 
    return *fMcPtrColl;
  }

  int DisplayBackgroundHits() { return fDisplayBackgroundHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void  SetKalRepPtrColl(const mu2e::KalRepPtrCollection** Coll) {
    fKalRepPtrColl = Coll;
  }

  void  SetListOfTracks(TObjArray* List) { fListOfTracks = List; }
  void  SetTrackColor  (Color_t      c ) { fTrackColor   = c   ; }

  void SetStrawHitColl(const mu2e::StrawHitCollection** Coll) { 
    fStrawHitColl = Coll;
  }

  void SetStrawHitPosColl (const mu2e::StrawHitPositionCollection** Coll) { 
    fStrawHitPosColl = Coll;
  }

  void SetStrawHitFlagColl(const mu2e::StrawHitFlagCollection** Coll) { 
    fStrawHitFlagColl = Coll;
  }

  void SetMcPtrColl(const mu2e::PtrStepPointMCVectorCollection** Coll) { 
    fMcPtrColl = Coll;
  }

  void SetStrawDigiMCColl(const mu2e::StrawDigiMCCollection** Coll) { 
    fStrawDigiMCColl = Coll;
  }

  void SetCalTimePeakColl(const mu2e::CalTimePeakCollection** Coll) { 
    fCalTimePeakColl = Coll;
  }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }

  void  SetDisplayBackgroundHits(Int_t Mode) { fDisplayBackgroundHits = Mode; }

//-----------------------------------------------------------------------------
// overloaded methods of TVisNode
//-----------------------------------------------------------------------------
  virtual int   InitEvent();
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* option = "");
          void  PaintXY (Option_t* option = "");
          void  PaintRZ (Option_t* option = "");
          void  PaintCal(Option_t* option = "");

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Clear(const char* Opt = "")       ; // **MENU**
  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TTrkVisNode,0)
};


#endif
