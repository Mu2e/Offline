///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdStrawTracker_hh
#define TEvdStrawTracker_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#ifndef __CINT__
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#else
namespace mu2e {
  class StrawHitCollection;
  class StrawHitPositionCollection;
  class StrawHitFlagCollection;
  class PtrStepPointMCVectorCollection;
  class TTracker;
};
#endif

class TEvdStation;

class TEvdStrawTracker: public TObject {
public:
  
protected:

  int        fNStations;
  TObjArray* fListOfStations;

  const mu2e::TTracker*  fTracker;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdStrawTracker(const mu2e::TTracker* Tracker = NULL);
  //  TEvdStrawTracker(const char* Name); 

  virtual ~TEvdStrawTracker();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int          NStations     () { return fNStations;      }
  TObjArray*   ListOfStations() { return fListOfStations; }

  TEvdStation* Station  (int I) { 
    return (TEvdStation*) fListOfStations->UncheckedAt(I); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* option = "");

  void  PaintXY   (Option_t* option = "");
  void  PaintRZ   (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdStrawTracker,0)
};


#endif
