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

#include "Stntuple/gui/TEvdCluster.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TEvdCluster)

//_____________________________________________________________________________
TEvdCluster::TEvdCluster(const mu2e::CaloCluster* Cl): TObject() {

  const int module_color[2] = {kRed, kMagenta};

  double xl, yl, rtrk, rcal;
  int    color;

  fCluster  = Cl;

  xl      = Cl->cog3Vector().x()+3904.1;
  yl      = Cl->cog3Vector().y();
  
  rtrk    = 50.*fCluster->energyDep()/105.; // radius od 

  fTrkEllipse = new TEllipse(xl,yl,rtrk,rtrk,0,360);
  fTrkEllipse->SetFillColor(2);
  fTrkEllipse->SetFillStyle(3001);

  art::ServiceHandle<mu2e::GeometryService> geom;

  if( geom->hasElement<mu2e::VaneCalorimeter>() ) {
    color = 2;
  } else if(geom->hasElement<mu2e::DiskCalorimeter>()){
    color = module_color[Cl->vaneId()];
  }

  fTrkEllipse->SetFillColor(color);

  rcal = 50;
  fCalEllipse = new TEllipse(xl,yl,rcal,rcal,0,360);
  fCalEllipse->SetFillColor(2);
  fCalEllipse->SetFillStyle(3001);

  fListOfCrystals = new TObjArray();
}

//-----------------------------------------------------------------------------
TEvdCluster::~TEvdCluster() {
  delete fTrkEllipse;
  delete fCalEllipse;
  delete fListOfCrystals;
}

//-----------------------------------------------------------------------------
void TEvdCluster::Paint(Option_t* Option) {
  const char oname[] = "TEvdCluster::Paint";

  int   iv;

  const char* view = TVisManager::Instance()->GetCurrentView();

  if      (strstr(view,"trkxy" ) != 0) {
    PaintXY(Option);
  }
  else if (strstr(view,"cal"   ) != 0) {
    sscanf(view,"cal,%i",&iv);
    if (iv == fCluster->vaneId()) {
      PaintCal(Option);
    }
  }
  else {
    printf("[%s] >>> ERROR: unknown view: %s, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TEvdCluster::PaintXY(Option_t* Option) {
  fTrkEllipse->Paint(Option);
}

//-----------------------------------------------------------------------------
// 2015-01-27: doesn't do anything
//-----------------------------------------------------------------------------
void TEvdCluster::PaintCal(Option_t* Option) {
  TEvdCrystal* cr; 
  THexagon     hex;
  
  int nc = fListOfCrystals->GetEntries();
  
  for (int i=0; i<nc; i++) {
    cr = Crystal(i);
    hex = *cr->Hexagon();
//     hex.SetFillStyle(1024);
//     hex.SetFillColor(0);
//     hex.SetLineColor(kBlue);
    //    hex.SetLineWidth(2);
    //    hex.Paint();
  }
  //  fCalEllipse->Paint(Option);
}


//_____________________________________________________________________________
Int_t TEvdCluster::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdCluster::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdCluster::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
void  TEvdCluster::Clear(Option_t* Option) {
  fListOfCrystals->Clear();
}

//_____________________________________________________________________________
void  TEvdCluster::Print(Option_t* Option) const {
  printf(" >>> WARNING: TEvdCluster::Print is not implemented yet\n");
}

