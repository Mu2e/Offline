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

// #include "art/Framework/Principal/Handle.h"

// #include "GeometryService/inc/GeometryService.hh"
// #include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/obj/TStnCrystal.hh"
#include "Stntuple/base/TVisManager.hh"
#include "Stntuple/base/THexagon.hh"
#include "Stntuple/base/TObjHandle.hh"
#include "Stntuple/obj/TCalHitData.hh"

// #include "CalorimeterGeom/inc/VaneCalorimeter.hh"
// #include "CalorimeterGeom/inc/Crystal.hh"
// #include "CalorimeterGeom/inc/Disk.hh"
// #include "CalorimeterGeom/inc/DiskCalorimeter.hh"
// #include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TStnCrystal)

//-----------------------------------------------------------------------------
TStnCrystal::TStnCrystal(THexIndex* Index, 
			 double     X0,
			 double     Y0,
			 double     Z0,
			 double     HexSize): TObject() {

  fHexIndex.fL = Index->fL;
  fHexIndex.fK = Index->fK;

  fCenter.SetXYZ(X0,Y0,Z0);
  fHexSize = HexSize;
  fHexagon.SetPos(X0,Y0);
  fHexagon.SetSize(30.);
  fHexagon.SetLineColor(1);
  fHexagon.SetFillColor(0);
  fHexagon.SetFillStyle(0);

  //  fDisk       = Disk;
  fListOfHits = new TObjArray(10);
  fNHits      = 0;
}

//-----------------------------------------------------------------------------
TStnCrystal::~TStnCrystal() {
  fListOfHits->Clear();
  delete fListOfHits;
}

//-----------------------------------------------------------------------------
void TStnCrystal::Paint(Option_t* Option) {
  const char oname[] = "TStnCrystal::Paint";

  //  int   iv;

  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) {
    PaintXY(Option);
  }
  if      (strstr(view,"cal"   ) != 0) {
//     sscanf(view,"cal,%i",&iv);
//     if (iv == fSectionNumber) {
    PaintCal(Option);
    //    }
  }
  else {
    printf("[%s] >>> ERROR: unknown view: %s, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
void TStnCrystal::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
void TStnCrystal::PaintCal(Option_t* Option) {
  fHexagon.Paint("f");
  fHexagon.Paint("same");
}


//_____________________________________________________________________________
Int_t TStnCrystal::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TStnCrystal::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TStnCrystal::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//-----------------------------------------------------------------------------
void TStnCrystal::Clear(Option_t* Opt) {
  SetFillStyle(0);
  SetFillColor(1);
  SetLineColor(1);
  fNHits     = 0;
  fEnergy    = 0.;

  fListOfHits->Clear();
}


//-----------------------------------------------------------------------------
void TStnCrystal::Print(Option_t* Opt) const {
  printf("TStnCrystal::Print is not implemented yet\n");
}
