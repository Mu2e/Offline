///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTrack_hh
#define TEvdTrack_hh

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
#else
namespace mu2e {
  class StrawHitCollection;
  class StrawHitPositionCollection;
  class StrawHitFlagCollection;
  class PtrStepPointMCVectorCollection;
};
#endif

class TEvdTrack: public TObject {
public:
  
protected:

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTrack();
  TEvdTrack(const char* Name); 

  virtual ~TEvdTrack();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintRZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  //  virtual void  Draw    (Option_t* option = "");
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");
  virtual void  Print(Option_t* Opt = "") const ; // **MENU**

  ClassDef(TEvdTrack,0)
};


#endif
