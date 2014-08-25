//_____________________________________________________________________________
// Feb 10 2001 P.Murat
//_____________________________________________________________________________
#include "TVirtualX.h"
#include "TPad.h"

#include "Stntuple/gui/THeaderVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

ClassImp(THeaderVisNode)

//_____________________________________________________________________________
THeaderVisNode::THeaderVisNode(const char* name, TStnHeaderBlock* h): 
  TVisNode(name) 
{
  fHeader = h;
  fText = new TText(0.1,0.1,"elki-palki");
  fText->SetTextFont(22);
}

//_____________________________________________________________________________
THeaderVisNode::~THeaderVisNode() {
  delete fText;
}


//_____________________________________________________________________________
int THeaderVisNode::InitEvent() {
  //
  TStnVisManager* vm = TStnVisManager::Instance();

  fHeader->fEventNumber = vm->Event()->event();
  fHeader->fRunNumber   = vm->Event()->run();

  return 0;
}

//_____________________________________________________________________________
void THeaderVisNode::Paint(Option_t* option) {
  //
				// parse option list
  PaintXY(option);
}


//_____________________________________________________________________________
void THeaderVisNode::PaintXY(Option_t* option) {
  // draw event/run

  char  text[100];
  sprintf(text,"Run = %7i  Event = %7i",
	  fHeader->RunNumber(), 
	  fHeader->EventNumber());

  fText->SetTextSize(0.4);
  fText->SetText(0.1,0.3,text);

  fText->Paint(option);
  gPad->Modified();
}


//_____________________________________________________________________________
void THeaderVisNode::PaintRZ(Option_t* option) {
}

//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t THeaderVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

