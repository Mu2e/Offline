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
#include "TBox.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TTrkVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

ClassImp(TTrkVisNode)

//_____________________________________________________________________________
TTrkVisNode::TTrkVisNode(const char* name, TStnTrackBlock* TrackBlock): TVisNode(name) {
  fTrackBlock = TrackBlock;
}

//_____________________________________________________________________________
TTrkVisNode::~TTrkVisNode() {
}


//_____________________________________________________________________________
void TTrkVisNode::Paint(Option_t* option) {
  //
				// parse option list

  if      (strstr(option,"XY" ) != 0) PaintXY (option);
  if      (strstr(option,"cal") != 0) PaintXY (option);
  else {
				// what is the default?
    printf(Form("[TTrkVisNode::Paint] >>> ERROR: Unknown option %s",option));
  }
}

//-----------------------------------------------------------------------------
int TTrkVisNode::InitEvent() {
  return 0;
}

//_____________________________________________________________________________
void TTrkVisNode::PaintXY(Option_t* option) {
  // draw tracks
 
     
//   double  d0, r, phi0, x0, y0;

//   //  TAnaDump::Instance()->printKalRep(0,"banner");

//   for (int i=0; i<fNTracks[0]; i++ ) {
//     trk = &(*fDem)[i];

//     d0   = trk->helix(0.).d0();
//     r    = fabs(1./trk->helix(0.).omega());
//     phi0 = trk->helix(0.).phi0();
    
//     x0   = -(r+d0)*sin(phi0);
//     y0   =  (r+d0)*cos(phi0);
	
//     e = new TEllipse(x0,y0,r);
//     e->SetFillStyle(3001);		// make it transparent

//     e->SetLineColor(kBlue-7);
//     list_of_ellipses.Add(e);
//     e->Draw();
//   }

  gPad->Modified();
}

//_____________________________________________________________________________
void TTrkVisNode::PaintRZ(Option_t* option) {
  // draw tracks
 
  gPad->Modified();
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
  //  static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TTrkVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

