///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TCanvas.h"
#include "TPad.h"
#include "TText.h"

//_____________________________________________________________________________
TCanvas* new_slide(const char* name, const char* slide_title, int nx, int ny,
		   int XSize, int YSize) {

  TCanvas* slide;

  if (name ) slide = new TCanvas(name,name,0,0,XSize,YSize);
  else       slide = new TCanvas("slide","slide",0,0,XSize,YSize);

  TPad *p1 = new TPad("p1", "p1",0.01,0.01,0.99,0.96);
  p1->Divide(nx,ny);
  p1->Draw();
  p1->Range(0,0,1,1);
  slide->cd();

  TPad *p2 = new TPad("p2", "p2",0.7,0.965,0.99,0.995);
  p2->Draw();
  p2->cd();

  time_t t2 = time(0);
  tm* t22 = localtime(&t2);
  TText *text = new TText(0.05,0.3,asctime(t22));
  text->SetTextFont(42);
  text->SetTextSize(0.5);
  text->Draw();
  p2->Modified();

  slide->cd();
  TPad *title = new TPad("title", "title",0.05,0.965,0.69,0.995);
  title->Draw();
  title->cd();

  if (slide_title) {
    TText *text = new TText(0.05,0.3,slide_title);
    text->SetTextFont(42);
    text->SetTextSize(0.5);
    text->Draw();
  }
  title->Modified();

  return slide;
}

