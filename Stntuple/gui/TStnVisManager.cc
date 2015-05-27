///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TVirtualX.h"

#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGStatusBar.h"
#include "TGaxis.h"
#include "TText.h"
#include "TGDoubleSlider.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGTextBuffer.h"
#include "TGLabel.h"

#include "Stntuple/gui/TEvdMainFrame.hh"

#include "Stntuple/gui/TTrkXYView.hh"
#include "Stntuple/gui/TTrkRZView.hh"
#include "Stntuple/gui/TCalView.hh"
#include "Stntuple/gui/TCrvView.hh"

#include "Stntuple/gui/TStnFrame.hh"
#include "Stntuple/gui/TStnVisManager.hh"


ClassImp(TStnVisManager)

//_____________________________________________________________________________
TStnVisManager::TStnVisManager(const char* Name, const char* Title) :
TVisManager(Name, Title)
{
	if (gROOT->IsBatch()) return;

	fMain = new  TEvdMainFrame(gClient->GetRoot(), 200, 100,
		kMainFrame | kVerticalFrame);
	//-----------------------------------------------------------------------------
	//  create menu bar
	//-----------------------------------------------------------------------------
	
	fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
	fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
	fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

	fMenu = new TGPopupMenu(gClient->GetRoot());
	fMenu->AddEntry("&Exit", M_EXIT);
	fMenu->Associate(fMain);

	fMenuHelp = new TGPopupMenu(gClient->GetRoot());
	fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
	fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
	fMenuHelp->AddSeparator();
	fMenuHelp->AddEntry("&About", M_HELP_ABOUT);
	fMenuHelp->Associate(fMain);

	fMenuBar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
	fMenuBar->AddPopup("&Menu", fMenu, fMenuBarItemLayout);
	fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

	fMain->AddFrame(fMenuBar, fMenuBarLayout);

	trkrBtnXY = new TGTextButton(fMain, "Tracker XY", kXYView);
	trkrBtnXY->Connect("Clicked()", "TStnVisManager", this, "HandleButtons()");
	trkrBtnXY->SetTextJustify(36);
	trkrBtnXY->SetMargins(0, 0, 0, 0);
	trkrBtnXY->SetWrapLength(-1);
	trkrBtnXY->MoveResize(16, 26, 98, 24);

	trkrBtnRZ = new TGTextButton(fMain, "Tracker RZ", kRZView);
	trkrBtnRZ->Connect("Clicked()", "TStnVisManager", this, "HandleButtons()");
	trkrBtnRZ->SetTextJustify(36);
	trkrBtnRZ->SetMargins(0, 0, 0, 0);
	trkrBtnRZ->SetWrapLength(-1);
	trkrBtnRZ->MoveResize(16, 58, 98, 24);

	calBtn = new TGTextButton(fMain, "Calorimeter", kCalView);
	calBtn->Connect("Clicked()", "TStnVisManager", this, "HandleButtons()");
	calBtn->SetTextJustify(36);
	calBtn->SetMargins(0, 0, 0, 0);
	calBtn->SetWrapLength(-1);
	calBtn->MoveResize(16, 90, 98, 24);

	crvBtn = new TGTextButton(fMain, "CRV", kCrvView);
	crvBtn->Connect("Clicked()", "TStnVisManager", this, "HandleButtons()");
	crvBtn->SetTextJustify(36);
	crvBtn->SetMargins(0, 0, 0, 0);
	crvBtn->SetWrapLength(-1);
	crvBtn->MoveResize(16, 122, 98, 24);

	timeWindowSlider = new TGDoubleHSlider(fMain, 100, kDoubleScaleBoth, TIMESLIDER_ID);
	timeWindowSlider->SetRange(0, 1695);
	timeWindowSlider->SetPosition(400, 1695);
	timeWindowSlider->MoveResize(150, 45, 200, 20);
	timeWindowSlider->Connect("PositionChanged()", "TStnVisManager", this, "HandleSlider()");
	
	timeWindowLowDisp = new TGTextEntry(fMain, timeWindowLowBuff = new TGTextBuffer(10), TIMELOW_DISP);
	timeWindowLowBuff->AddText(0, "400");
	timeWindowLowDisp->MoveResize(150, 70, 40, 20);
	timeWindowLowDisp->Connect("ReturnPressed()", "TStnVisManager", this, "HandleText()");

	timeWindowHighDisp = new TGTextEntry(fMain, timeWindowHighBuff = new TGTextBuffer(10), TIMEHIGH_DISP);
	timeWindowHighBuff->AddText(0, "1695");
	timeWindowHighDisp->MoveResize(310, 70, 40, 20);
	timeWindowHighDisp->Connect("ReturnPressed()", "TStnVisManager", this, "HandleText()");

	TGLabel *sliderLabelLow = new TGLabel(fMain, "0");
	sliderLabelLow->SetTextJustify(36);
	sliderLabelLow->SetMargins(0, 0, 0, 0);
	sliderLabelLow->SetWrapLength(-1);
	sliderLabelLow->MoveResize(140, 25, 30, 20);

	TGLabel *sliderLabelHigh = new TGLabel(fMain, "1695");
	sliderLabelHigh->SetTextJustify(36);
	sliderLabelHigh->SetMargins(0, 0, 0, 0);
	sliderLabelHigh->SetWrapLength(-1);
	sliderLabelHigh->MoveResize(330, 25, 30, 20);

	updaterBtn = new TGTextButton(fMain, "Update", UPDATER_BTN);
	updaterBtn->Connect("Clicked()", "TStnVisManager", this, "HandleButtons()");
	updaterBtn->SetTextJustify(36);
	updaterBtn->SetMargins(0, 0, 0, 0);
	updaterBtn->SetWrapLength(-1);
	updaterBtn->MoveResize(220, 120, 60, 20);

	//-----------------------------------------------------------------------------
	// views
	//-----------------------------------------------------------------------------
	fTrkXYView = new TTrkXYView();
	fTrkRZView = new TTrkRZView();
	// for vanes - 4 "views"
	fCalView[0] = new TCalView(0);
	fCalView[1] = new TCalView(1);

	fCrvView[0] = new TCrvView(0); //right
	fCrvView[0]->SetTimeWindow(0, 1695);
	fCrvView[1] = new TCrvView(1); //left
	fCrvView[1]->SetTimeWindow(0, 1695);
	fCrvView[2] = new TCrvView(2); //topds
	fCrvView[2]->SetTimeWindow(0, 1695);
	fCrvView[3] = new TCrvView(3); //downstream
	fCrvView[3]->SetTimeWindow(0, 1695);
	fCrvView[4] = new TCrvView(4); //upstream
	fCrvView[4]->SetTimeWindow(0, 1695);
	fCrvView[5] = new TCrvView(8); //topts
	fCrvView[5]->SetTimeWindow(0, 1695);
	

	fListOfDetectors = new TObjArray(10);

	//-----------------------------------------------------------------------------
	// final actions
	//-----------------------------------------------------------------------------
	fMain->MapSubwindows();
	fMain->Resize(fMain->GetDefaultSize());
	fMain->Resize(400, 150);

	fMain->SetWindowName(Title);

	fMain->MapWindow();

	fMinStation = 0;
	fMaxStation = 50;
	// by default, no timing constraints
	fTMin = 0;
	fTMax = 1.e5;

	fTimePeak = -1;
}

//_____________________________________________________________________________
TStnVisManager::~TStnVisManager() {

	if (!gROOT->IsBatch()) {

		// delete tracking views
		delete fTrkXYView;
		delete fTrkRZView;
		// only two views for disk calorimeter
		delete fCalView[0];
		delete fCalView[1];

		delete fCrvView[0];
		delete fCrvView[1];
		delete fCrvView[2];
		delete fCrvView[3];
		delete fCrvView[4];
		delete fCrvView[5];

		delete fMenuBarHelpLayout;

		delete fMenuBarItemLayout;
		delete fMenu;

		delete fMenuBarLayout;
		delete fMenuBar;

		delete fMain;

		delete fListOfDetectors;
	}
}

//_____________________________________________________________________________
TStnVisManager* TStnVisManager::Instance() {
	if (TVisManager::fgInstance != NULL) {
		return (TStnVisManager*) TVisManager::fgInstance;
	}
	else {
		return new TStnVisManager();
	}
}

//_____________________________________________________________________________
TCanvas* TStnVisManager::NewCanvas(const char* Name,
	const char* Title,
	Int_t       SizeX,
	Int_t       SizeY)
{
	TStnFrame* win = new TStnFrame(Name, Title, 0, SizeX, SizeY);
	TCanvas*c = win->GetCanvas();
	DeclareCanvas(c);
	return c;
}


//_____________________________________________________________________________
Int_t TStnVisManager::OpenTrkXYView() {
	// open new XY view of the detector with the default options

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "xy_view_%i", n);
	sprintf(title, "XY view number %i", n);

	TStnFrame* win = new TStnFrame(name, title, kXYView, 740, 760);
	TCanvas* c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	p1->Range(-800., -800., 800., 800.);
	p1->cd();
	fTrkXYView->Draw();

	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenTrkXYView(TTrkXYView* mother, Axis_t x1, Axis_t y1,
	Axis_t x2, Axis_t y2)
{
	// open new XY view of the detector with the default options

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "xy_view_%i", n);
	sprintf(title, "XY view number %i", n);

	// try to preserve the aspect ratio
	Int_t   xsize, ysize;

	xsize = 540;
	ysize = (Int_t) (xsize*TMath::Abs((y2 - y1) / (x2 - x1)) + 20);

	TStnFrame* win = new TStnFrame(name, title, kXYView, xsize, ysize);
	TCanvas* c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	p1->Range(x1, y1, x2, y2);
	p1->cd();
	mother->Draw();

	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}

//-----------------------------------------------------------------------------
// open new RZ view of the detector with the default options
//-----------------------------------------------------------------------------
Int_t TStnVisManager::OpenTrkRZView() {

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "rz_view_%i", n);
	sprintf(title, "RZ view number %i", n);

	TStnFrame* win = new TStnFrame(name, title, kRZView, 1500, 500);
	TCanvas* c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	//  p1->Range(8500.,-200.,12500.,800.);
	p1->Range(-2000., -300., 3000., 700.);
	p1->cd();
	fTrkRZView->Draw();

	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}

//-----------------------------------------------------------------------------
// open new RZ view of the detector with the default options
//-----------------------------------------------------------------------------
Int_t TStnVisManager::OpenTrkRZView(TTrkRZView* mother,
	Axis_t      x1, Axis_t y1,
	Axis_t      x2, Axis_t y2)
{

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "rz_view_%i", n);
	sprintf(title, "RZ view number %i", n);
	//-----------------------------------------------------------------------------
	// try to preserve the aspect ratio
	//-----------------------------------------------------------------------------
	Int_t   xsize, ysize;

	xsize = 540;
	ysize = (Int_t) (xsize*TMath::Abs((y2 - y1) / (x2 - x1)) + 20);

	TStnFrame* win = new TStnFrame(name, title, kRZView, xsize, ysize);
	TCanvas* c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	p1->Range(x1, y1, x2, y2);
	p1->cd();
	mother->Draw();

	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}


//_____________________________________________________________________________
Int_t TStnVisManager::OpenCalView() {
	// open new calorimeter view of the detector with the default options
	// start from the disk calorimeter

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "cal_view_%i", n);
	sprintf(title, "CAL view number %i", n);

	TStnFrame* win = new TStnFrame(name, title, kCalView, 1200, 600);
	TCanvas*   c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	//-----------------------------------------------------------------------------
	// the disk calorimeter view has two pads, one per disk
	// the vane-based calorimeter display should have 4 pads in this view
	//-----------------------------------------------------------------------------
	// divide horisontally
	p1->Divide(2, 1);
	// ranges in mm
	p1->cd(1);
	gPad->Range(-800., -800., 800., 800.);
	fCalView[0]->SetPad(gPad);
	fCalView[0]->Draw("cal");
	gPad->Modified();

	p1->cd(2);
	gPad->Range(-800., -800., 800., 800.);
	fCalView[1]->SetPad(gPad);
	fCalView[1]->Draw("cal");
	gPad->Modified();
	// draw title
	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenCalView(TObject* Mother, Axis_t x1, Axis_t y1,
	Axis_t x2, Axis_t y2)
{
	// open new calorimeter view of the detector with the default options

	//   int n = fListOfCanvases->GetSize();

	//   char name[100], title[100];

	//   sprintf(name,"ces_view_%i",n);
	//   sprintf(title,"CES view number %i",n);

	// 				// try to preserve the aspect ratio
	//   Int_t   xsize, ysize;

	//   xsize = 540;
	//   ysize = (Int_t) (xsize*TMath::Abs((y2-y1)/(x2-x1))+20);

	//   TStnFrame* win = new TStnFrame(name, title, kCesStripView, xsize,ysize);
	//   TCanvas* c = win->GetCanvas();
	//   fListOfCanvases->Add(c);
	//   c->Divide(1,2);

	//   TString name1(name);
	//   name1 += "_1";
	//   TPad* p1 = (TPad*) c->FindObject(name1);

	//   p1->Divide(2,1);

	//   p1->cd(1);
	//   gPad->Range(x1,y1,x2,y2);

	//   fCalSectionView[0]->Draw();

	//   p1->cd(2);
	//   gPad->Range(x1,y1,x2,y2);
	//   fCalSectionView[1]->Draw();

	//   TString name_title(name);
	//   name1 += "_title";
	//   TPad* title_pad = (TPad*) c->FindObject(name_title);
	//   title_pad->cd();
	//   fTitleNode->Draw();

	//   c->Modified();
	//   c->Update();
	return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenCrvView() {
	TText Tl;
	Tl.SetTextSize(0.2);

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "crv_view_%i", n);
	sprintf(title, "CRV view number %i", n);

	TStnFrame* win = new TStnFrame(name, title, kCrvView, 2000, 600);
	TCanvas*   c = win->GetCanvas();
	c->SetFixedAspectRatio(kTRUE);
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);

	p1->Divide(1, 10, 0.003, 0.003);
	// ranges in mm
	p1->cd(1);
	gPad->Range(2698., -6570., 18750., -6415.); // Right CRV sans TS region
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[0]->SetPad(gPad);
	fCrvView[0]->Draw("crv");
	Tl.DrawText(3000, -6565, "RIGHT");
	gPad->SetEditable(kFALSE);
	gPad->Modified();
		
	p1->cd(2);
	gPad->Range(2698., -1420., 18750., -1260.); // Left CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[1]->SetPad(gPad);
	fCrvView[1]->Draw("crv");
	Tl.DrawText(3000, -1420, "LEFT");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(3);
	gPad->Range(2698., 2560., 18750., 2715.); // Top DS CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[2]->SetPad(gPad);
	fCrvView[2]->Draw("crv");
	Tl.DrawText(3000, 2560, "TOP DS");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//Axis Pad
	p1->cd(4);
	gPad->Range(2698., 0., 18750., 10.);
	TGaxis *a1 = new TGaxis(2698., 7., 18750., 7., 2698., 18750., 50510, "");
	a1->SetName("Zaxis");
	a1->SetTitle("Z (mm)");
	a1->SetLabelSize(0.3);
	a1->SetTitleSize(0.3);
	a1->SetLabelOffset(0.1);
	a1->SetTitleOffset(0.3);
	a1->SetTickSize(0.1);
	a1->Draw();
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(5);
	gPad->Range(-2220., -6570., 2750., -6415.); // Right CRV TS region
	gPad->SetFixedAspectRatio(kTRUE);
	//fCrvView[0]->SetPad(gPad);
	fCrvView[0]->Draw("crv");
	Tl.DrawText(-2150, -6565, "RIGHT");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(6);
	gPad->Range(-2220., 2560., 2750., 2715.); // Top TS CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[5]->SetPad(gPad);
	fCrvView[5]->Draw("crv");
	Tl.DrawText(-2150, 2560, "TOP TS");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//Axis Pad
	p1->cd(7);
	gPad->Range(-2220., 0., 2750., 10.);
	TGaxis *a2 = new TGaxis(-2220., 7., 2750., 7., -2220., 2750., 50510, "");
	a2->SetName("Zaxis2");
	a2->SetTitle("Z (mm)");
	a2->SetLabelSize(0.3);
	a2->SetTitleSize(0.3);
	a2->SetLabelOffset(0.1);
	a2->SetTitleOffset(0.3);
	a2->SetTickSize(0.1);
	a2->Draw();
	gPad->SetEditable(kFALSE);
	gPad->Modified();	
		

	p1->cd(8);
	gPad->Range(-2220., 18695., 2780., 18970.); // Downstream CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[3]->SetPad(gPad);
	fCrvView[3]->Draw("crv");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(9);
	gPad->Range(-2220., -2415., 2780., -2260.); // Upstream CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[4]->SetPad(gPad);
	fCrvView[4]->Draw("crv");
	Tl.DrawText(-2000., -2300., "DWNSTRM");
	Tl.DrawText(50., -2340., "UPSTRM");
	gPad->SetEditable(kFALSE);
	gPad->Modified();


	p1->cd(10);
	gPad->Range(-2220., 0., 2800., 10.);
	TGaxis *a3 = new TGaxis(-2220., 7., 2780., 7., -2220, 2780, 50510, "");
	a3->SetName("Yaxis");
	a3->SetTitle("Y (mm)");
	a3->SetLabelSize(0.3);
	a3->SetTitleSize(0.3);
	a3->SetLabelOffset(0.2);
	a3->SetTitleOffset(0.3);
	a3->SetTickSize(0.1);
	a3->Draw();
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	// draw title
	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();
	title_pad->SetEditable(kFALSE);

	c->cd();
	gPad->SetEditable(kFALSE);
	c->Modified();
	c->Update();
	return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenCrvView(TCrvView* mother, Axis_t x1, Axis_t y1,
	Axis_t x2, Axis_t y2)
{

	int n = fListOfCanvases->GetSize();

	char name[100], title[100];

	sprintf(name, "crv_view_%i", n);
	sprintf(title, "CRV view number %i", n);

	// try to preserve the aspect ration
	Int_t   xsize, ysize;

	xsize = 540;
	ysize = (Int_t) (xsize*TMath::Abs((y2 - y1) / (x2 - x1)) + 20);

	TStnFrame* win = new TStnFrame(name, title, kCrvView, xsize, ysize);
	TCanvas* c = win->GetCanvas();
	fListOfCanvases->Add(c);

	TString name1(name);
	name1 += "_1";
	TPad* p1 = (TPad*) c->FindObject(name1);
	p1->Range(x1, y1, x2, y2);
	p1->cd();
	mother->Draw();

	TString name_title(name);
	name1 += "_title";
	TPad* title_pad = (TPad*) c->FindObject(name_title);
	title_pad->cd();
	fTitleNode->Draw();

	c->Modified();
	c->Update();
	return 0;
}

void TStnVisManager::UpdateViews()
{
	TIter it(fListOfCanvases);
	while (TCanvas* c = (TCanvas*) it.Next()) {
		TIter it1(c->GetListOfPrimitives());
		while (TObject* o = it1.Next()) {
			if (o->InheritsFrom("TPad")) {
				TPad* pad = (TPad*) o;
				MarkModified(pad);
			}
		}
		c->Modified();
		c->Update();
	}
}


//_____________________________________________________________________________
void TStnVisManager::CloseWindow() {
	// Called when window is closed via the window manager.

	delete this;
}


//-----------------------------------------------------------------------------
void TStnVisManager::SetStations(int IMin, int IMax) {
	fMinStation = IMin;
	fMaxStation = IMax;
}

//-----------------------------------------------------------------------------
void TStnVisManager::SetTimePeak(int I) {
	fTimePeak = I;
}

//_____________________________________________________________________________
void TStnVisManager::HandleButtons()
{
	// Handle different buttons.

	TGButton *btn = (TGButton *) gTQSender;
	int id = btn->WidgetId();

	switch (id)
	{
	case kXYView:
		OpenTrkXYView();
		break;
	case kRZView:
		OpenTrkRZView();
		break;
	case kCalView:
		OpenCalView();
		break;
	case kCrvView:
		OpenCrvView();
		break;
	case UPDATER_BTN:
		for (unsigned int i = 0; i < 6; i++)
		{
			fCrvView[i]->SetTimeWindow(timeWindowSlider->GetMinPosition(), timeWindowSlider->GetMaxPosition());
		}
		UpdateViews();
		break;
	default:
		printf("Unknown button clicked\n");
		break;
	}
}

//_____________________________________________________________________________
void TStnVisManager::HandleSlider()
{
	// Handle slider widget

	Int_t id;
	TGFrame *frm = (TGFrame *) gTQSender;
	TGDoubleSlider *sd = (TGDoubleSlider *) frm;
	id = sd->WidgetId();

	switch (id) {
		//case TStnVisManager::TIMESLIDER_ID:
		// Update text boxes with max and min values
		
	case TIMESLIDER_ID:
		timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
		gClient->NeedRedraw(timeWindowLowDisp);

		timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
		gClient->NeedRedraw(timeWindowHighDisp);
		break;
	default:
		break;
	}
}

//_____________________________________________________________________________
void TStnVisManager::HandleText() //const char * /*text*/)
{
	// Handle text entry widgets

	TGTextEntry *te = (TGTextEntry *) gTQSender;
	Int_t id = te->WidgetId();

	float textBoxNum;

	switch (id) {
	case TIMELOW_DISP:
		try{
			textBoxNum = boost::lexical_cast<float>(timeWindowLowDisp->GetText());
			if (textBoxNum < 0 || textBoxNum > timeWindowSlider->GetMaxPosition())
				timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
			else
			{
				timeWindowSlider->SetPosition(textBoxNum, timeWindowSlider->GetMaxPosition());
				timeWindowLowDisp->SetText(boost::lexical_cast<std::string>(textBoxNum).c_str());
			}
		}
		catch (boost::bad_lexical_cast &){
			timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
		}		
		break;
	case TIMEHIGH_DISP:
		try{
			textBoxNum = boost::lexical_cast<float>(timeWindowHighDisp->GetText());
			if (textBoxNum > 1695 || textBoxNum < timeWindowSlider->GetMinPosition())
				timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
			else
			{
				timeWindowSlider->SetPosition(timeWindowSlider->GetMinPosition(), textBoxNum);
				timeWindowHighDisp->SetText(boost::lexical_cast<std::string>(textBoxNum).c_str());
			}
		}
		catch (boost::bad_lexical_cast &){
			timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
		}
		break;
	default:
		break;
	}
}