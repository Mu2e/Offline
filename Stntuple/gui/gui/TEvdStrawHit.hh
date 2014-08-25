///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdStrawHit_hh
#define TEvdStrawHit_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TVector2.h"
#include "TLine.h"

#ifndef __CINT__

#include "RecoDataProducts/inc/StrawHit.hh"

#else

namespace mu2e {
  class StrawHit;
};

#endif

class TEvdStrawHit: public TObject {
public:

  enum { 
    kInTimeBit     = 0x1 << 0,
    kConversionBit = 0x1 << 1
  };

protected:

  const mu2e::StrawHit*  fHit;
  int        fMask;			// hit mask
  int        fColor;
  TVector2   fPos;
  TVector2   fStrawDir;
  double     fSigW;      		// error in the wire direction
  double     fSigR;      		// error in radial direction
  TLine      fLineW;
  TLine      fLineR;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdStrawHit() {}
  TEvdStrawHit(const mu2e::StrawHit* Hit);

  virtual ~TEvdStrawHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::StrawHit*  StrawHit() { return fHit; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetMask (int Mask ) { fMask = Mask ;}
  void SetColor(int Color) { 
    fLineW.SetLineColor(Color); 
    fLineR.SetLineColor(Color); 
  }
  void SetSigR(double Sig) { fSigR = Sig; }
  void SetSigW(double Sig) { fSigW = Sig; }

  void SetPos(double X, double Y) { fPos.Set(X,Y); }
  void SetStrawDir(double X, double Y) { fStrawDir.Set(X,Y); }

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint   (Option_t* option = "");
  virtual void  PaintXY   (Option_t* option = "");
  virtual void  PaintCal   (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdStrawHit,0)
};


#endif
