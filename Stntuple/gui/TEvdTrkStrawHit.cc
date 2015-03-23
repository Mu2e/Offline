///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TMath.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TObjArray.h"

#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TEvdTrkStrawHit.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(TEvdTrkStrawHit)

//_____________________________________________________________________________
TEvdTrkStrawHit::TEvdTrkStrawHit(const mu2e::TrkStrawHit* Hit): TObject() {
  fHit = Hit;

  fLineW.SetX1(fPos.X()-fStrawDir.X()*fSigW);
  fLineW.SetY1(fPos.Y()-fStrawDir.Y()*fSigW);
  fLineW.SetX2(fPos.X()+fStrawDir.X()*fSigW);
  fLineW.SetY2(fPos.Y()+fStrawDir.Y()*fSigW);

  fLineR.SetX1(fPos.X()+fStrawDir.Y()*fSigR);
  fLineR.SetY1(fPos.Y()-fStrawDir.X()*fSigR);
  fLineR.SetX2(fPos.X()-fStrawDir.Y()*fSigR);
  fLineR.SetY2(fPos.Y()+fStrawDir.X()*fSigR);

  const mu2e::Straw* straw = &Hit->straw();
  
  double zw     = straw->getMidPoint().z();
  double rw     = straw->getMidPoint().perp();
  double rdrift = fHit->driftRadius();

  fEllipse.SetX1(zw);
  fEllipse.SetY1(rw);
  fEllipse.SetR1(rdrift);
  fEllipse.SetR2(rdrift);
  fEllipse.SetFillStyle(3001);
  fEllipse.SetFillColor(kBlue+2);
}

//-----------------------------------------------------------------------------
TEvdTrkStrawHit::~TEvdTrkStrawHit() {
}

//-----------------------------------------------------------------------------
void TEvdTrkStrawHit::Paint(Option_t* Option) {
  const char oname[] = "TEvdTrkStrawHit::Paint";

  const char* view = TVisManager::Instance()->GetCurrentView();
//-----------------------------------------------------------------------------
// define lines
//-----------------------------------------------------------------------------
  if      (strstr(view,"trkxy" ) != 0) PaintXY (Option);
  if      (strstr(view,"trkrz" ) != 0) PaintRZ (Option);
  else if (strstr(view,"cal"   ) != 0) PaintCal(Option);
  else {
    printf("[%s] >>> ERROR: unknown view: %s, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TEvdTrkStrawHit::PaintXY(Option_t* Option) {
  fLineW.Paint(Option);
  fLineR.Paint(Option);
}

//_____________________________________________________________________________
void TEvdTrkStrawHit::PaintRZ(Option_t* Option) {

  if (fHit->isActive()) fEllipse.SetFillColor(kRed-9 );
  else                  fEllipse.SetFillColor(kBlue+2);

  fEllipse.Paint(Option);
}

//_____________________________________________________________________________
void TEvdTrkStrawHit::PaintCal(Option_t* option) {
  // nothing to draw...
}


//_____________________________________________________________________________
Int_t TEvdTrkStrawHit::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdTrkStrawHit::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdTrkStrawHit::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

