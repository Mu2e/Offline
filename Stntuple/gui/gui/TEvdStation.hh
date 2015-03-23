///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdStation_hh
#define TEvdStation_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#ifndef __CINT__
#include "TrackerGeom/inc/Device.hh"
#else
namespace mu2e {
  class   Device;
};
#endif

class TEvdStation;
class TEvdPanel;

class TEvdStation: public TObject {
public:
  
protected:
  int                   fID;
  int                   fNPanels;
  TObjArray*            fListOfPanels;

  const mu2e::Device*  fDevice;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdStation();
  TEvdStation(int ID, const mu2e::Device* Device); 

  virtual ~TEvdStation();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int          NPanels     () { return fNPanels;      }
  TObjArray*   ListOfPanels() { return fListOfPanels; }

  TEvdPanel* Panel  (int I) { 
    return (TEvdPanel*) fListOfPanels->UncheckedAt(I); 
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

  ClassDef(TEvdStation,0)
};


#endif
