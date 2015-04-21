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
#include "TVector3.h"
#include "TLine.h"
#include "TEllipse.h"

#ifndef __CINT__

#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#else
namespace mu2e {
  class StrawHit;
  class StrawDigiMC;
  class PtrStepPointMCVector;
};
#endif

class TEvdStraw;

class TEvdStrawHit: public TObject {
public:

  enum { 
    kInTimeBit     = 0x1 << 0,
    kConversionBit = 0x1 << 1
  };

protected:
					// backward pointers to the reconstruction objects

  const mu2e::StrawHit*             fHit;
  const mu2e::StrawDigiMC*          fStrawDigiMC;
  const mu2e::PtrStepPointMCVector* fListOfStepPointMCs;

  TEvdStraw*  fStraw;                 // pointer to the straw -  geometry

  int         fMask;			// hit mask
  int         fColor;
  TVector3    fPos;			// position in 3D, Z=zwire
  TVector2    fDir;                     // direction of the straw
  double      fSigW;      		// error in the wire direction
  double      fSigR;      		// error in radial direction
  TLine       fLineW;			// paint on XY view
  TLine       fLineR;
  TEllipse    fEllipse;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdStrawHit();

  TEvdStrawHit(const mu2e::StrawHit*        Hit,
	       TEvdStraw*                   Straw,
	       const mu2e::StrawDigiMC*     StrawDigiMC,
	       double X, double Y, double Z, 
	       double                       Wx,
	       double                       Wy,
	       double                       SigW,
	       double                       SigR,
	       int                          Mask, 
	       int                          Color);

  virtual ~TEvdStrawHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::StrawHit*        StrawHit()    { return fHit; }
  const mu2e::StrawDigiMC*     StrawDigiMC() { return fStrawDigiMC; }
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

  void SetPos(double X, double Y, double Z) { fPos.SetXYZ(X,Y,Z); }
  void SetStrawDir(double X, double Y) { fDir.Set(X,Y); }

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

  ClassDef(TEvdStrawHit,0)
};


#endif
