///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTrkStrawHit_hh
#define TEvdTrkStrawHit_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TVector2.h"
#include "TLine.h"

#ifndef __CINT__
#include "KalmanTests/inc/TrkStrawHit.hh"
#else
namespace mu2e {
  class TrkStrawHit;
};
#endif

class TEvdStraw;

class TEvdTrkStrawHit: public TObject {
public:

  enum { 
    kInTimeBit     = 0x1 << 0,
    kConversionBit = 0x1 << 1
  };

protected:

  const mu2e::TrkStrawHit*  fHit;
  TEvdStraw*                fStraw;

  int        fMask;			// hit mask
  int        fColor;
  TVector2   fPos;
  TVector2   fStrawDir;
  double     fSigW;      		// error in the wire direction
  double     fSigR;      		// error in radial direction
  TLine      fLineW;
  TLine      fLineR;
  TEllipse   fEllipse;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTrkStrawHit() {}
  TEvdTrkStrawHit(const mu2e::TrkStrawHit* Hit);

  virtual ~TEvdTrkStrawHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::TrkStrawHit*  TrkStrawHit() { return fHit; }
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

  virtual void  Paint      (Option_t* option = "");
  virtual void  PaintXY    (Option_t* option = "");
  virtual void  PaintRZ    (Option_t* option = "");
  virtual void  PaintCal   (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TEvdTrkStrawHit,0)
};


#endif
