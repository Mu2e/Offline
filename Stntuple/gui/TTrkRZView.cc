///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TObjArray.h"

#include "Stntuple/base/TVisNode.hh"

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

  vm->SetCurrentView("trkxy");

  fCenter->Paint(Option);

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node =  vm->GetNode(i);
    node->PaintXY(Option);
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
  Axis_t x, y, x1, x2, y1, y2;

  switch (event) {

  case kButton1Down:

    printf(" TTrkRZView::ExecuteEvent  kButton1Down\n");

    fPx1 = px;
    fPy1 = py;
    fPx2 = px;
    fPy2 = py;

    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, TVirtualX::kHollow);

    break;

  case kButton1Motion:
				// redraw the opaque rectangle
    gPad->SetCursor(kCross);

    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, TVirtualX::kHollow);
    fPx2 = px;
    fPy2 = py;
    gVirtualX->DrawBox(fPx1, fPy1, fPx2, fPy2, TVirtualX::kHollow);

    break;
  case kButton2Up:
    printf(" TTrkRZView::ExecuteEvent kButton2Up\n");
    break;
  case kButton2Down:
    printf("TTrkRZView::ExecuteEvent  kButton2Down\n");
    break;
  case kButton3Up:
    printf(" TTrkRZView::ExecuteEvent kButton3Up\n");
    break;
  case kButton3Down:
    printf(" TTrkRZView::ExecuteEvent kButton3Down\n");
    break;
  case kButton1Up:
    printf(" TTrkRZView::ExecuteEvent kButton1Up\n");

				        // open new window only if something
					// has really been selected (it is a 
					// rectangle, but not a occasional 
					// click)

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
  case 26:
    gPad->GetRange(x1,y1,x2,y2);
    gPad->Range(x1*1.1,y1*1.1,x2*1.1,y2*1.1);
    gPad->Modified();
    break;
  case 60:
    gPad->GetRange(x1,y1,x2,y2);
    gPad->Range(x1/1.1,y1/1.1,x2/1.1,y2/1.1);
    gPad->Modified();
    break;
  default:
    //    printf(" ----- event = %i\n",event);
    break;
  }
}

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

