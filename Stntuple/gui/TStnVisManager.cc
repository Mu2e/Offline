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
#include "TGButton.h"


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

	fMenuSubdetectors = new TGPopupMenu(gClient->GetRoot());
	fMenuSubdetectors->AddEntry("&Tracker(XY view)", M_TRACKER_XY);
	fMenuSubdetectors->AddEntry("&Tracker(RZ view)", M_TRACKER_RZ);
	fMenuSubdetectors->AddEntry("&Calorimeter", M_CALORIMETER_XY);
	fMenuSubdetectors->AddEntry("&CRV", M_CRV_XY);
	fMenuSubdetectors->AddEntry("&Exit", M_EXIT);
	fMenuSubdetectors->Associate(fMain);

	fMenuHelp = new TGPopupMenu(gClient->GetRoot());
	fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
	fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
	fMenuHelp->AddSeparator();
	fMenuHelp->AddEntry("&About", M_HELP_ABOUT);
	fMenuHelp->Associate(fMain);

	fMenuBar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
	fMenuBar->AddPopup("&Subdetectors", fMenuSubdetectors, fMenuBarItemLayout);
	fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

	fMain->AddFrame(fMenuBar, fMenuBarLayout);

	TGTextButton *trkrBtn = new TGTextButton(fMain, "Tracker", 1);
	trkrBtn->Connect("Clicked()", "TEvdMainFrame", this, "HandleButtons()");
	trkrBtn->SetTextJustify(36);
	trkrBtn->SetMargins(0, 0, 0, 0);
	trkrBtn->SetWrapLength(-1);
	trkrBtn->Resize(98, 24);
	fMain->AddFrame(trkrBtn, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
	trkrBtn->MoveResize(16, 16, 98, 24);

	TGTextButton *calBtn = new TGTextButton(fMain, "Calorimeter", 2);
	calBtn->Connect("Clicked()", "TEvdMainFrame", this, "HandleButtons()");
	calBtn->SetTextJustify(36);
	calBtn->SetMargins(0, 0, 0, 0);
	calBtn->SetWrapLength(-1);
	calBtn->Resize(98, 24);
	fMain->AddFrame(calBtn, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
	calBtn->MoveResize(16, 48, 98, 24);

	TGTextButton *crvBtn = new TGTextButton(fMain, "CRV", 3);
	crvBtn->Connect("Clicked()", "TEvdMainFrame", this, "HandleButtons()");
	crvBtn->SetTextJustify(36);
	crvBtn->SetMargins(0, 0, 0, 0);
	crvBtn->SetWrapLength(-1);
	crvBtn->Resize(98, 24);
	fMain->AddFrame(crvBtn, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
	crvBtn->MoveResize(16, 80, 98, 24);

	//TGDoubleHSlider *timeWindowSlider = new TGDoubleHSlider(fMain, 100, kDoubleScaleBoth, TIMESLIDER_ID);
	//timeWindowSlider->SetRange(0, 1695);
	////timeWindowSlider->Connect("PositionChanged()", "TEvdMainFrame", this, "DoSlider()");

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
	fMain->Resize(200, 100);

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
		delete fMenuSubdetectors;

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

	p1->Divide(1, 8, 0.005, 0.005);
	// ranges in mm
	p1->cd(1);
	gPad->Range(-2500., -6570., 19000., -6415.); // Right CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[0]->SetPad(gPad);
	fCrvView[0]->Draw("crv");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//p1->cd(2);
	//gPad->Range(2500.,-1395.,19000.,-1240.); // Left CRV
	//gPad->SetFixedAspectRatio(kTRUE);
	//fCrvView[1]->SetPad(gPad);
	//fCrvView[1]->Draw("crv");
	//gPad->SetEditable(kFALSE);
	//gPad->Modified();

	p1->cd(2);
	gPad->Range(-2500., -1395., 19000., -1240.); // Left CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[1]->SetPad(gPad);
	fCrvView[1]->Draw("crv");
	Tl.SetTextSize(0.2);
	Tl.DrawText(0, -1290, "RIGHT");
	Tl.DrawText(2000, -1320, "LEFT");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//p1->cd(3);
	//gPad->Range(2350.,2585.,19000.,2740.); // Top DS CRV
	//gPad->SetFixedAspectRatio(kTRUE);
	//fCrvView[2]->SetPad(gPad);
	//fCrvView[2]->Draw("crv");
	//gPad->SetEditable(kFALSE);
	//gPad->Modified(); 

	p1->cd(3);
	gPad->Range(-2500., 2585., 19000., 2740.); // Top DS CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[2]->SetPad(gPad);
	fCrvView[2]->Draw("crv");
	Tl.DrawText(2000, 2665, "TOP DS");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(4); //p1->cd(6);
	gPad->Range(-2500., 2585., 19000., 2740.); // Top TS CRV
	//gPad->Range(-2200.,2585.,2800.,2740.);
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[5]->SetPad(gPad);
	fCrvView[5]->Draw("crv");
	Tl.DrawText(3000, 2665, "TOP TS");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//Axis Pad
	p1->cd(5);
	gPad->Range(-2500., 0., 19000., 10.);
	TGaxis *a1 = new TGaxis(-2500., 7., 19000., 7., -2500, 19000, 50510, "");
	a1->SetName("Zaxis");
	a1->SetTitle("Z (mm)");
	a1->SetLabelSize(0.3);
	a1->SetTitleSize(0.3);
	a1->SetLabelOffset(0.1);
	a1->SetTitleOffset(0.3);
	a1->SetTickSize(0.1);
	a1->Draw();
	p1->cd(5);
	p1->cd(5);
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(6); //p1->cd(4);
	gPad->Range(200., 18755., 2800., 18910.); // Downstream CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[3]->SetPad(gPad);
	fCrvView[3]->Draw("crv");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(7); //p1->cd(5);
	gPad->Range(200., -2415., 2800., -2260.); // Upstream CRV
	gPad->SetFixedAspectRatio(kTRUE);
	fCrvView[4]->SetPad(gPad);
	fCrvView[4]->Draw("crv");
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	p1->cd(8);
	gPad->Range(200., 0., 2800., 10.);
	TGaxis *a2 = new TGaxis(200., 7., 2800., 7., 200, 2800, 50510, "");
	a2->SetName("Yaxis");
	a2->SetTitle("Y (mm)");
	a2->SetLabelSize(0.3);
	a2->SetTitleSize(0.3);
	a2->SetLabelOffset(0.2);
	a2->SetTitleOffset(0.3);
	a2->SetTickSize(0.1);
	a2->Draw();
	gPad->SetEditable(kFALSE);
	gPad->Modified();

	//p1->cd(7); //p1->cd(6);
	//gPad->Range(-2200.,2585.,2800.,2740.); // Top TS CRV
	//gPad->SetFixedAspectRatio(kTRUE);
	//fCrvView[5]->SetPad(gPad);
	//fCrvView[5]->Draw("crv");
	//gPad->SetEditable(kFALSE);
	//gPad->Modified();

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

