#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TArrow.h"
#include "TCut.h"
#include "TBox.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "Math/Math.h"
#include "THStack.h"

class HelixDiag  {
  public:
    HelixDiag(TTree* hdiag) : _hdiag(hdiag) {}

    void CenterRes();
    void CenterPos();

    TTree* _hdiag;
};

void HelixDiag::CenterRes() {
  TH2F* crcomp = new TH2F("crcomp","Reco vs true Center Radius;Reco Center Radius (mm);MC Center Radius(mm)",50,150.0,450.0,50,150.0,450.0);
  TH1F* crres = new TH1F("crres","Center Radius Resolution;#Delta R (mm)",100,-150.0,150.0);
  TH2F* cfcomp = new TH2F("cfcomp","Reco vs true Center #Phi;Reco Center #phi;MC Center #phi",50,-3.15,3.15,50,-3.15,3.15);
  TH1F* cfres = new TH1F("cfres","Center Radius #Phi Resolution;#Delta #phi",100,-0.4,0.4);

  _hdiag->Project("crcomp","sqrt(mch._center.dx^2+mch._center.dy^2):sqrt(rhel._center.dx^2+rhel._center.dy^2)");
  _hdiag->Project("crres","sqrt(mch._center.dx^2+mch._center.dy^2)-sqrt(rhel._center.dx^2+rhel._center.dy^2)");
  _hdiag->Project("cfcomp","atan2(mch._center.dy,mch._center.dx):atan2(rhel._center.dy,rhel._center.dx)");
  _hdiag->Project("cfres","atan2(mch._center.dy,mch._center.dx)-atan2(rhel._center.dy,rhel._center.dx)");

  TCanvas* crcan = new TCanvas("crcan","Center Resolution",800,800);
  crcan->Divide(2,2);
  crcan->cd(1);
  crcomp->Draw();
  crcan->cd(2);
  crres->Draw();
  crcan->cd(3);
  cfcomp->Draw();
  crcan->cd(4);
  cfres->Draw();
}

void HelixDiag::CenterPos() {
  TH2F* rcpos = new TH2F("rcpos","Reco Center Position; x(mm) y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* mccpos = new TH2F("mccpos","MC Center Position; x(mm) y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* xcomp = new TH2F("xcomp","Reco vs true Center x;Reco Center x (mm);MC Center x (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* ycomp = new TH2F("ycomp","Reco vs true Center y;Reco Center y (mm);MC Center y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  _hdiag->Project("rcpos","rhel._center.dy:rhel._center.dx");
  _hdiag->Project("mccpos","mch._center.dy:mch._center.dx");
  _hdiag->Project("xcomp","rhel._center.dx:mch._center.dx");
  _hdiag->Project("ycomp","rhel._center.dy:mch._center.dy");

  TCanvas* cpcan = new TCanvas("cpcan","Center Position",800,800);
  cpcan->Divide(2,2);
  cpcan->cd(1);
  rcpos->Draw();
  cpcan->cd(2);
  mccpos->Draw();
  cpcan->cd(3);
  xcomp->Draw();
  cpcan->cd(4);
  ycomp->Draw();
}
