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

#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdFace.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TEvdStrawTracker.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"


ClassImp(TEvdStrawTracker)

//_____________________________________________________________________________
TEvdStrawTracker::TEvdStrawTracker(const mu2e::TTracker* Tracker): TObject() {

  const mu2e::Station* station;
  TEvdStation*         s;

  fTracker = Tracker;

  if (Tracker == NULL) {
    printf(">>> TEvdStrawTracker::TEvdStrawTracker ERROR: Tracker = NULL\n");
  }

  fNStations      = fTracker->nStations();
  fListOfStations = new TObjArray(fNStations);

  for (int i=0; i<fNStations; i++) {
    station = &fTracker->getStations().at(i);

    s = new TEvdStation(i,station);
    fListOfStations->Add(s);
  }
}

//-----------------------------------------------------------------------------
// need destructor...
//-----------------------------------------------------------------------------
TEvdStrawTracker::~TEvdStrawTracker() {
  if (fListOfStations) {
    fListOfStations->Delete();
    delete fListOfStations;
  }
}

//-----------------------------------------------------------------------------
void TEvdStrawTracker::Paint(Option_t* option) {

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
void TEvdStrawTracker::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdStrawTracker::PaintRZ(Option_t* option) {
  // draw tracker

  TEvdStation*   station;
  TEvdPlane*     plane;
  TEvdFace*      face;
  TEvdPanel*     panel;
  TEvdStraw*     straw;

  int            nplanes, nfaces, npanels, nlayers, nstraws;

  for (int ist=0; ist<fNStations; ist++) {
    station = Station(ist);
    nplanes = station->NPlanes();
    for (int ipl=0; ipl<nplanes; ipl++) {
      plane  = station->Plane(ipl);
      nfaces = plane->NFaces();
      for (int ifc=0; ifc<nfaces; ifc++) {
	face    = plane->Face(ifc);
	npanels = face->NPanels();
	for (int ipanel=0; ipanel<npanels; ipanel++) {
	  panel   = face->Panel(ipanel);
	  nlayers = panel->NLayers();
	  for (int il=0; il<nlayers; il++) {
	    nstraws = panel->NStraws(il);
	    for (int is=0; is<nstraws; is++) {
	      straw  = panel->Straw(il,is);
	      straw->PaintRZ();
	    }
	  }
	}
      }
    }
  }
}

//_____________________________________________________________________________
Int_t TEvdStrawTracker::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdStrawTracker::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdStrawTracker::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

