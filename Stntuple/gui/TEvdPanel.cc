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

#include "TrackerGeom/inc/Layer.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Sector.hh"

#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TStnVisManager.hh"


ClassImp(TEvdPanel)

//_____________________________________________________________________________
TEvdPanel::TEvdPanel(): TObject() {
}

//_____________________________________________________________________________
TEvdPanel::TEvdPanel(int ID, const mu2e::Sector* Sector, TEvdStation* Station): TObject() {

  TEvdStraw* evd_straw;

  fID      = ID;
  fNLayers = Sector->nLayers();
  fSector   = Sector;
					// assume that the number of straws is the same
  int id, id0;

  id0 = fID*fNLayers*2;
  for (int il=0; il<fNLayers; il++) {
    const mu2e::Layer* layer = &fSector->getLayer(il);

    fNStraws[il] = layer->nStraws();
    fListOfStraws[il] = new TObjArray(fNStraws[il]);

    for (int is=0; is<fNStraws[il]; is++) {
      const mu2e::Straw* straw = &layer->getStraw(is);
      id        = straw->index().asInt();
      evd_straw = new TEvdStraw(id,straw,this);
      fListOfStraws[il]->Add(evd_straw);
    }
    id0 +=fNStraws[il];
  }

}

//_____________________________________________________________________________
TEvdPanel::~TEvdPanel() {
}

//-----------------------------------------------------------------------------
void TEvdPanel::Paint(Option_t* option) {


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
void TEvdPanel::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdPanel::PaintRZ(Option_t* option) {
  // draw straws
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

