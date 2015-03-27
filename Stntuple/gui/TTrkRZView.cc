///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TObjArray.h"

#include "Stntuple/base/TVisNode.hh"

#include "Stntuple/gui/TStnFrame.hh"
#include "Stntuple/gui/TTrkRZView.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(TTrkRZView)

//_____________________________________________________________________________
TTrkRZView::TTrkRZView() {
  fCenter = new TMarker(0.,0,kPlus);
  fCenter->SetMarkerColor(kBlue);
  fCenter->SetMarkerSize(3.);
}

//_____________________________________________________________________________
TTrkRZView::~TTrkRZView() {
  delete fCenter;
}


//_____________________________________________________________________________
void TTrkRZView::Paint(Option_t* Option) {
  //
  TStnVisManager* vm = TStnVisManager::Instance();

  vm->SetCurrentView("trkrz");

  fCenter->Paint(Option);

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node =  vm->GetNode(i);
    node->Paint(Option);
  }
  gPad->Modified();
}

//_____________________________________________________________________________
Int_t TTrkRZView::DistancetoPrimitive(Int_t px, Int_t py) {
  //

  Int_t min_dist = 9999;
  Int_t dist;

  TStnVisManager* vm = TStnVisManager::Instance();

  vm->SetClosestObject(NULL,9999);
  vm->SetClosestDetElement(NULL,9999);
//-----------------------------------------------------------------------------
// find closest object and the corresponding distance
//-----------------------------------------------------------------------------
//  TDetectorElement::SetClosest(NULL,9999);

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node = vm->GetNode(i);
    dist = node->DistancetoPrimitiveXY(px,py);
    if (dist < min_dist) {
      min_dist = dist;
				// closest object may be one of the managed by
				// the node

      vm->SetClosestObject(node->GetClosestObject(),dist);
    }
  }
//-----------------------------------------------------------------------------
// find closest detector element (not sure if I really need this)
//-----------------------------------------------------------------------------
//  if (TDetectorElement::GetMinDist() < vm->GetMinDistDetElement()) {
//     vm->SetClosestDetElement(TDetectorElement::GetClosest(),
// 			     TDetectorElement::GetMinDist());
//  }

  if (vm->GetMinDist() > 5) vm->SetClosestObject(this,0);
//-----------------------------------------------------------------------------
// prepare output
//-----------------------------------------------------------------------------
  gPad->SetSelected(vm->GetClosestObject());
  return vm->GetMinDist();
}


//_____________________________________________________________________________
void TTrkRZView::ExecuteEvent(Int_t event, Int_t px, Int_t py) {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  TObject* o = vm->GetClosestObject();
  if (o && (o != this) && (vm->GetMinDist() < 5) ) {
    o->ExecuteEvent(event,px,py);
    return;
  }
//-----------------------------------------------------------------------------
// view...
//-----------------------------------------------------------------------------
  Axis_t x, y, x1, x2, y1, y2, dx, dy;

  double     xx,yy;
  //  int        px, py;
  
  double shift_scale = 0.02;
  double zoom_scale  = 0.05;

  TVirtualX::EBoxMode ebox_mode = TVirtualX::kHollow;
  //  TVirtualX::EBoxMode ebox_mode = TVirtualX::kFilled;

  gVirtualX->SetLineColor(1);
  gVirtualX->SetFillColor(1);
  gVirtualX->SetFillStyle(3003);

  switch (event) {

  case kButton1Down:

    if (vm->DebugLevel() > 0) printf(" TTrkRZView::ExecuteEvent  kButton1Down\n");

    fPx1 = px;
    fPy1 = py;
    fPx2 = px;
    fPy2 = py;

    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, ebox_mode);

    break;

  case kButton1Motion:
				// redraw the opaque rectangle
    gPad->SetCursor(kCross);

    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, ebox_mode);
    fPx2 = px;
    fPy2 = py;
    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, ebox_mode);

    break;
  case kKeyPress:
    if (vm->DebugLevel() > 0) printf(" TTrkRZView::ExecuteEvent kKeyPress: px=%3i py:%i\n",px,py);

    if (px == py) {
      Axis_t  x1new, x2new, y1new, y2new;

      gPad->GetRange(x1,y1,x2,y2);

      if (char(px) == 'z') {            // zoom in
	x1new = x1+(x2-x1)/2.*zoom_scale;
	x2new = x2-(x2-x1)/2.*zoom_scale;
	y1new = y1+(y2-y1)/2.*zoom_scale;
	y2new = y2-(y2-y1)/2.*zoom_scale;
      }
      else if (char(px) == 'Z') { 	// zoom out
	x1new = x1-(x2-x1)/2.*zoom_scale;
	x2new = x2+(x2-x1)/2.*zoom_scale;
	y1new = y1-(y2-y1)/2.*zoom_scale;
	y2new = y2+(y2-y1)/2.*zoom_scale;
      }

      gPad->Range(x1new,y1new,x2new,y2new);
      gPad->Modified();
      gPad->Update();
    }
    break;
  case kButton2Up:
    if (vm->DebugLevel() > 0) printf(" TTrkRZView::ExecuteEvent kButton2Up\n");
    break;
  case kButton2Down:
    if (vm->DebugLevel() > 0) printf("TTrkRZView::ExecuteEvent  kButton2Down\n");
    break;
  case kButton3Up:
    if (vm->DebugLevel() > 0) printf(" TTrkRZView::ExecuteEvent kButton3Up\n");
    break;
  case kButton3Down:
    if (vm->DebugLevel() > 0) printf(" TTrkRZView::ExecuteEvent kButton3Down\n");
    break;
  case kButton1Up:
    printf(" TTrkRZView::ExecuteEvent kButton1Up\n");
//-----------------------------------------------------------------------------
// open new window only if something has really been selected (it is a 
// rectangle, but not a occasional click)
//-----------------------------------------------------------------------------
    if ( (fPx1 != fPx2) && (fPy1 != fPy2) ) {
      x1 = gPad->AbsPixeltoX(fPx1);
      y1 = gPad->AbsPixeltoY(fPy1);
      x2 = gPad->AbsPixeltoX(fPx2);
      y2 = gPad->AbsPixeltoY(fPy2);

      if (x1 > x2) {
	x  = x1;
	x1 = x2;
	x2 = x;
      }

      if (y1 > y2) {
	y  = y1;
	y1 = y2;
	y2 = y;
      }
      vm->OpenTrkRZView(this,x1,y1,x2,y2);
    }
    break;
  case kArrowKeyPress: // 25
//-----------------------------------------------------------------------------
// arrow key released
//-----------------------------------------------------------------------------
    if (vm->DebugLevel() > 0) printf(" ARROW key pressed  = %3i  py = %3i\n",px,py);
    fCursorX = px;
    fCursorY = py;
//     gPad->GetRange(x1,y1,x2,y2);
//     gPad->Range(x1*1.1,y1*1.1,x2*1.1,y2*1.1);
//    gPad->Modified();
    break;
  case kArrowKeyRelease: // 26
//-----------------------------------------------------------------------------
// arrow key released
//-----------------------------------------------------------------------------
    if (vm->DebugLevel() > 0) printf(" ARROW key  releasedx = %3i  py = %3i\n",px,py);

    gPad->GetRange(x1,y1,x2,y2);

    if (px == fCursorX) {
					         // up or down arrow
      dy = y2-y1;
      if (fCursorY > py) {
					         // down arrow

	gPad->Range(x1,y1-shift_scale*dy,x2,y2-shift_scale*dy);
      }
      else {
	gPad->Range(x1,y1+shift_scale*dy,x2,y2+shift_scale*dy);
      }
    }
    else if (py == fCursorY) {
					// left or right arrow
      dx = x2-x1;
      if (fCursorX > px) {
					         // left arrow
	gPad->Range(x1+shift_scale*dx,y1,x2+shift_scale*dx,y2);
      }
      else {
	gPad->Range(x1-shift_scale*dx,y1,x2-shift_scale*dx,y2);
      }
    }
    gPad->Modified();
    gPad->Update();

    break;
  default:
    if (vm->DebugLevel() > 0) printf(" ----- event = %i px:%4i py:%4i\n",
				     event,px,py);
    break;
  }

//   if (event != kMouseLeave) {                  // signal was already emitted for this event
//     DrawEventStatus(event, px, py);
//   }
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
// void TTrkRZView::DrawEventStatus(event,px,py) {

//   if (!TestBit(kShowEventStatus) || !selected) return;
 
//   TVirtualPad* savepad;
//   savepad = gPad;

//   gPad = GetSelectedPad();
 
//   fCanvasImp->SetStatusText(selected->GetTitle(),0);
//   fCanvasImp->SetStatusText(selected->GetName(),1);

//   snprintf(atext, kTMAX, "%d,%d", px, py);

//   fCanvasImp->SetStatusText(atext,2);
//   fCanvasImp->SetStatusText(selected->GetObjectInfo(px,py),3);
//}

//-----------------------------------------------------------------------------
void    TTrkRZView::SetStations(int I1, int I2) {
  TStnVisManager* vm = TStnVisManager::Instance();

  vm->SetStations(I1, I2);
}

//-----------------------------------------------------------------------------
void    TTrkRZView::SetTimePeak(int I) {
  TStnVisManager* vm = TStnVisManager::Instance();

  vm->SetTimePeak(I);
}


//-----------------------------------------------------------------------------
char* TTrkRZView::GetObjectInfo(Int_t Px, Int_t Py) const {
  // need to find the active frame:

  double xx, yy;

  static char info[200];

  xx = gPad->AbsPixeltoX(Px);
  yy = gPad->AbsPixeltoY(Py);
  
  sprintf(info,"Z=%9.3f R=%8.3f",xx,yy);

  return info;
  
}
