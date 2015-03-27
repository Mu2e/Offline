///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdStraw_hh
#define TEvdStraw_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

namespace mu2e {
  class Straw;
};

class TEvdStation;
class TEvdStraw;
class TEvdStrawHit;
class TEvdPanel;

class TEvdStraw: public TObject {
public:
  
protected:
  int                fIndex;
  int                fNHits;
  TObjArray*         fListOfHits;

  TEvdPanel*         fPanel; 		// backward pointer
  TArc*              fArc;              // circle ...

  const mu2e::Straw* fStraw;

  TObjArray*         fListOfStepPointMCs;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdStraw();
  TEvdStraw(int Index, const mu2e::Straw* Straw, TEvdPanel* Panel); 

  virtual ~TEvdStraw();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int NHits() { return fNHits; }

  TEvdStrawHit*   Hit(int I) { 
    return (TEvdStrawHit*) fListOfHits->UncheckedAt(I); 
  }

  const mu2e::Straw*   GetStraw() { return fStraw; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void AddHit(TObject* Hit) { fListOfHits->Add(Hit); }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY (Option_t* option = "");
  virtual void  PaintRZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  virtual void Clear(const char* Opt = "")       ;
  virtual void Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdStraw,0)
};


#endif
