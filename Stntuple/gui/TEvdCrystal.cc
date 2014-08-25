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

#include "Stntuple/gui/TEvdCrystal.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "Stntuple/base/THexagon.hh"
#include "Stntuple/base/TObjHandle.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TEvdCrystal)

//-----------------------------------------------------------------------------
  TEvdCrystal::TEvdCrystal(const mu2e::Crystal* Cr, TDisk* Disk): TObject() {
  fCrystal = Cr;

  const CLHEP::Hep3Vector  *pos;

  pos = &Cr->position();

  fHexagon.SetPos(pos->x(),pos->y());
  fHexagon.SetSize(30.);
  fHexagon.SetLineColor(1);
  fHexagon.SetFillColor(0);
  fHexagon.SetFillStyle(0);

  fDisk       = Disk;
  fListOfHits = new TClonesArray("TObjHandle",10);
  fNHits      = 0;
}

//-----------------------------------------------------------------------------
TEvdCrystal::~TEvdCrystal() {
}

//-----------------------------------------------------------------------------
void TEvdCrystal::Paint(Option_t* Option) {
  const char oname[] = "TEvdCrystal::Paint";

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
void TEvdCrystal::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
void TEvdCrystal::PaintCal(Option_t* Option) {
  fHexagon.Paint("f");
  fHexagon.Paint("same");
}


//_____________________________________________________________________________
Int_t TEvdCrystal::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdCrystal::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdCrystal::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//-----------------------------------------------------------------------------
void TEvdCrystal::AddHit(const mu2e::CaloCrystalHit* CrystalHit) {
  fEnergy  += CrystalHit->energyDep();
  TObjHandle* h = new ((*fListOfHits)[fNHits]) TObjHandle((void*) CrystalHit);
  fNHits++;
}




//-----------------------------------------------------------------------------
void TEvdCrystal::Clear(Option_t* Opt) {
  SetFillStyle(0);
  SetFillColor(1);
  SetLineColor(1);
  fNHits     = 0;
  fEnergy    = 0.;

  fListOfHits->Clear();
}


//-----------------------------------------------------------------------------
void TEvdCrystal::Print(Option_t* Opt) const {

  TObjHandle*                  h;
  const mu2e::CaloCrystalHit*  hit;

  printf (" X0 = %10.3f Y0 = %10.3f  E = %10.3f  njits = %5i\n",
	  X0(),Y0(),fEnergy,fNHits);

  printf("----------------------------------------------------------------\n");
  printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
  printf("----------------------------------------------------------------\n");
  for (int i=0; i<fNHits; i++) {

    h   = (TObjHandle*) fListOfHits->At(i);
    hit = (const mu2e::CaloCrystalHit*) h->Object();

    printf("%7i  %10.3f %10.3f %10.3f %5i\n",
	   hit->id(),
	   hit->time(),
	   hit->energyDep(),
	   hit->energyDepTotal(),
	   hit->numberOfROIdsUsed());
  }
}
