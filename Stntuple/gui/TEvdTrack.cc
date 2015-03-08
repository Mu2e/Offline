///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
//
// BaBar interface:
// ----------------
//      r     = fabs(1./om);
//      phi0  = Trk->helix(0.).phi0();
//      x0    =  -(1/om+d0)*sin(phi0);
//      y0    =   (1/om+d0)*cos(phi0);
//
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

#include "Stntuple/gui/TEvdTrack.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "Stntuple/base/TObjHandle.hh"


ClassImp(TEvdTrack)

//-----------------------------------------------------------------------------
TEvdTrack::TEvdTrack(): TObject() {
}

//-----------------------------------------------------------------------------
TEvdTrack::~TEvdTrack() {
}

//-----------------------------------------------------------------------------
void TEvdTrack::Paint(Option_t* Option) {
  const char oname[] = "TEvdTrack::Paint";

  //  int   iv;

  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) {
    PaintXY(Option);
  }
  else if (strstr(view,"trkrz" ) != 0) {
    PaintRZ(Option);
  }
  else if (strstr(view,"cal"   ) != 0) {
//-----------------------------------------------------------------------------
// calorimeter-specific view: do not draw tracks
//-----------------------------------------------------------------------------
  }
  else {
    printf("[%s] >>> ERROR: unknown view: %s, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
void TEvdTrack::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
void TEvdTrack::PaintRZ(Option_t* Option) {
}

//_____________________________________________________________________________
Int_t TEvdTrack::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdTrack::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdTrack::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}


//-----------------------------------------------------------------------------
void TEvdTrack::Clear(Option_t* Opt) {
  //  SetLineColor(1);
}


//-----------------------------------------------------------------------------
void TEvdTrack::Print(Option_t* Opt) const {

//   TObjHandle*                  h;
//   const mu2e::CaloCrystalHit*  hit;

//   printf (" X0 = %10.3f Y0 = %10.3f  E = %10.3f  njits = %5i\n",
// 	  X0(),Y0(),fEnergy,fNHits);

//   printf("----------------------------------------------------------------\n");
//   printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
//   printf("----------------------------------------------------------------\n");
//   for (int i=0; i<fNHits; i++) {

//     h   = (TObjHandle*) fListOfHits->At(i);
//     hit = (const mu2e::CaloCrystalHit*) h->Object();

//     printf("%7i  %10.3f %10.3f %10.3f %5i\n",
// 	   hit->id(),
// 	   hit->time(),
// 	   hit->energyDep(),
// 	   hit->energyDepTotal(),
// 	   hit->numberOfROIdsUsed());
//   }
}
