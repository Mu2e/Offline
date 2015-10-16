///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdFace_hh
#define TEvdFace_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

namespace mu2e {
  class Face;
};

class TEvdPlane;
class TEvdPanel;

class TEvdFace: public TObject {
public:
  
protected:
  int                 fID;
  int                 fNPanels;         // a face has 3 panels in it
  TObjArray*          fListOfPanels;

  const TEvdPlane*    fPlane; 		// backward pointers
  const mu2e::Face*   fFace;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdFace();
  TEvdFace(int Number, const mu2e::Face* Face, const TEvdPlane* Plane); 

  virtual ~TEvdFace();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int        NPanels     () { return fNPanels;  }
  TEvdPanel* Panel  (int I) { return (TEvdPanel*) fListOfPanels->UncheckedAt(I); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint   (Option_t* option = "");
          void  PaintXY (Option_t* option = "");
          void  PaintRZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdFace,0)
};


#endif
