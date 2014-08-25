///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"

#include "Stntuple/base/TVisNode.hh"
//#include "TGeant/TSubdetector.hh"

#include "Stntuple/gui/TCalView.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TCalVisNode.hh"

ClassImp(TCalView)

//_____________________________________________________________________________
TCalView::TCalView(int Section): TNamed("CalView","CalView") {
  fSectionToDisplay = Section;
  fPad              = 0;
}

//_____________________________________________________________________________
TCalView::~TCalView() {
}


//_____________________________________________________________________________
void TCalView::Paint(Option_t* Option) {
// cal view has two pads how does it know about it?

  TString view("cal");

  TString opt = Option;

  TStnVisManager* vm = TStnVisManager::Instance();

  view += Form(",%i",fSectionToDisplay);

  vm->SetCurrentView(view);

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node =  vm->GetNode(i);
    node->Paint(opt);
  }

  int np = gPad->GetListOfPrimitives()->GetEntries();
  for (int i=0; i<np; i++) {
    TObject* o = gPad->GetListOfPrimitives()->At(i);
    if (o->InheritsFrom("TPad")) {
      ((TPad*) o)->Modified();
    }
  }

  gPad->Modified();
  gPad->Update();
}

//_____________________________________________________________________________
Int_t TCalView::DistancetoPrimitive(Int_t px, Int_t py) {
  //

  Int_t min_dist = 9999;
  Int_t dist;

  TStnVisManager* vm = TStnVisManager::Instance();

  vm->SetCurrentView(Form("cal,%i",fSectionToDisplay));

  vm->SetClosestObject(NULL,9999);
  vm->SetClosestDetElement(NULL,9999);
//-----------------------------------------------------------------------------
// find closest object and the corresponding distance
//-----------------------------------------------------------------------------
//  TDetectorElement::SetClosest(NULL,9999);

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node = vm->GetNode(i);
    dist = node->DistancetoPrimitive(px,py);
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
//   if (TDetectorElement::GetMinDist() < vm->GetMinDistDetElement()) {
//     vm->SetClosestDetElement(TDetectorElement::GetClosest(),
// 			     TDetectorElement::GetMinDist());
//   }

  if (vm->GetMinDist() > 5) vm->SetClosestObject(this,0);
//-----------------------------------------------------------------------------
// prepare output
//-----------------------------------------------------------------------------
  gPad->SetSelected(vm->GetClosestObject());
  return vm->GetMinDist();
}


//_____________________________________________________________________________
void TCalView::ExecuteEvent(Int_t event, Int_t px, Int_t py) {
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

    printf(" TCalView::ExecuteEvent  kButton1Down\n");

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
    printf(" kButton2Up\n");
    break;
  case kButton2Down:
    printf("kButton2Down\n");
    break;
  case kButton3Up:
    printf(" kButton3Up\n");
    break;
  case kButton3Down:
    printf(" kButton3Down\n");
    break;
  case kButton1Up:
    printf(" TCalView::ExecuteEvent kButton1Up\n");

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
      vm->OpenCalView(this,x1,y1,x2,y2);
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


//_____________________________________________________________________________
void TCalView::SetMinClusterEnergy(float MinEnergy) {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node = vm->GetNode(i);

    if (node->InheritsFrom("TCalVisNode")) {
      TCalVisNode* cvn = (TCalVisNode*) node;
      cvn->SetMinClusterEnergy(MinEnergy);
      cvn->InitEvent();
    }
  }
}

//_____________________________________________________________________________
void TCalView::SetMinCrystalEnergy(float MinEnergy) {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  Int_t n = vm->GetNNodes();
  for (int i=0; i<n; i++) {
    TVisNode* node = vm->GetNode(i);

    if (node->InheritsFrom("TCalVisNode")) {
      TCalVisNode* cvn = (TCalVisNode*) node;
      cvn->SetMinCrystalEnergy(MinEnergy);
      cvn->InitEvent();
    }
  }
}


//_____________________________________________________________________________
void TCalView::PrintClosestCrystal() {
  //

  TStnVisManager* vm = TStnVisManager::Instance();

  TObject* clo = vm->GetClosestObject();

  if (clo && clo->InheritsFrom("TEvdCrystal")) {
    clo->Print();
  }
  else {
    printf("TCalView::PrintClosestCrystal ERROR: closest object is not a crystal\n");
  }

}
