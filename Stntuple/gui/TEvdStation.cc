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

#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TStnVisManager.hh"


ClassImp(TEvdStation)

//_____________________________________________________________________________
TEvdStation::TEvdStation(): TObject() {
  fListOfPanels = NULL;

}

//_____________________________________________________________________________
TEvdStation::TEvdStation(int ID, const mu2e::Device* Device): TObject() {

  int         id;
  TEvdPanel*  panel;

  fID      = ID;
  fDevice  = Device;
  fNPanels = Device->nSectors();

  fListOfPanels = new TObjArray(fNPanels);

  for (int i=0; i<fNPanels; i++) {
    const mu2e::Sector* sector = &fDevice->getSector(i);

    id    = fNPanels*Device->id()+i;
    panel = new TEvdPanel(id,sector,this);

    fListOfPanels->Add(panel);
  }
}

//_____________________________________________________________________________
TEvdStation::~TEvdStation() {
  fListOfPanels->Delete();
  delete fListOfPanels;
}

//-----------------------------------------------------------------------------
void TEvdStation::Paint(Option_t* option) {
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
void TEvdStation::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdStation::PaintRZ(Option_t* option) {
  // draw calorimeter
}

//_____________________________________________________________________________
Int_t TEvdStation::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdStation::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdStation::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

