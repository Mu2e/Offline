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
#include "TObjArray.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdFace.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "TTrackerGeom/inc/Face.hh"


ClassImp(TEvdFace)

//_____________________________________________________________________________
TEvdFace::TEvdFace(): TObject() {
  fID           = -1;
  fListOfPanels = NULL;
  fNPanels      = 0;
}

//_____________________________________________________________________________
TEvdFace::TEvdFace(int ID, const mu2e::Face* Face, const TEvdPlane* Plane): TObject() {

  int         id;
  TEvdPanel*  evd_panel;

  fID      = ID;
  fPlane   = Plane;
  fFace    = Face;
  fNPanels = Face->nPanels();

  fListOfPanels = new TObjArray(fNPanels);

  for (int i=0; i<fNPanels; i++) {
    const mu2e::Panel* panel = &fFace->getPanels().at(i);

    id        = -1; // fNPanels*Face->id()+i;
    evd_panel = new TEvdPanel(id,panel,this);

    fListOfPanels->Add(evd_panel);
  }
}

//_____________________________________________________________________________
TEvdFace::~TEvdFace() {
  fListOfPanels->Delete();
  delete fListOfPanels;
}

//-----------------------------------------------------------------------------
void TEvdFace::Paint(Option_t* option) {
  // paints one disk (.. or vane, in the past), i.e. section

				// parse option list

  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) PaintXY (option);
  else if (strstr(view,"trkrz" ) != 0) PaintRZ (option);
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  gPad->Modified();
}


//_____________________________________________________________________________
void TEvdFace::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdFace::PaintRZ(Option_t* option) {
}

//_____________________________________________________________________________
Int_t TEvdFace::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdFace::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdFace::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

