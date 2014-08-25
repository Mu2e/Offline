///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TStrawHitVisNode_hh
#define TStrawHitVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__
#include "RecoDataProducts/inc/StrawHitCollection.hh"
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
  class CalTimePeakCollection;
  class CalTimePeak;
};
#endif

class TStrawHitVisNode: public TVisNode {
public:
  enum {
    kPickHits     = 0,
    kPickTracks   = 1,
    kPickClusters = 2
  };
  
protected:
  //  TObjArray**    fListOfClusters;

  const mu2e::StrawHitCollection**             fStrawHitColl;
  const mu2e::StrawHitPositionCollection**     fStrawHitPosColl;  //
  const mu2e::StrawHitFlagCollection**         fStrawHitFlagColl; //
  const mu2e::CalTimePeakCollection**          fCalTimePeakColl;  //
  const mu2e::PtrStepPointMCVectorCollection** fMcPtrColl; 
 
  TArc*         fArc;

  TObjArray*    fListOfStrawHits;

  const mu2e::CalTimePeak*  fTimePeak;

  Int_t         fDisplayBackgroundHits;
  Int_t         fTimeWindow;
  Int_t         fPickMode;
  Int_t         fUseStereoHits;
  double        fEventTime;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TStrawHitVisNode() {}
  TStrawHitVisNode(const char* Name); 

  virtual ~TStrawHitVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
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

  void SetCalTimePeakColl(const mu2e::CalTimePeakCollection** Coll) { 
    fCalTimePeakColl = Coll;
  }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }

  void  SetDisplayBackgroundHits(Int_t Mode) { fDisplayBackgroundHits = Mode; }

  //  virtual void  Draw    (Option_t* option = "");

  int InitEvent();

  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintRZ (Option_t* option = "");
  virtual void  PaintCal(Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TStrawHitVisNode,0)
};


#endif
