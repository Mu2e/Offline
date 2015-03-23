///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdPanel_hh
#define TEvdPanel_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

namespace mu2e {
  class Sector;
};

class TEvdStation;
class TEvdStraw;

class TEvdPanel: public TObject {
public:
  
protected:
  int        fID;
  int        fNLayers;
  int        fNStraws[2];
  TObjArray* fListOfStraws[2];

  TEvdStation*         fStation; 		// backward pointer
  const mu2e::Sector*  fSector;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdPanel();
  TEvdPanel(int Number, const mu2e::Sector* Sector, TEvdStation* Station); 

  virtual ~TEvdPanel();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int          NLayers     () { return fNLayers;  }
  int          NStraws     (int I) { return fNStraws[I];  }

  TEvdStraw* Straw  (int Layer, int I) { 
    return (TEvdStraw*) fListOfStraws[Layer]->UncheckedAt(I); 
  }
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

  ClassDef(TEvdPanel,0)
};


#endif
