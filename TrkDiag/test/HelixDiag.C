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
#include "TFitResult.h"
#include "Math/Math.h"
#include "THStack.h"
#include "TGraphErrors.h"

class HelixDiag  {
  public:
    HelixDiag(TTree* hdiag) : _hdiag(hdiag), _helixOK("helixOK"), _mchelixOK("mchelixOK"),
    _hhout("hh._outlier"), _stereo("hh._stereo"), _tdiv("hh._tdiv"),
    _tdivonly("hh._tdiv && !hh._stereo"), _resphi("hh._resphi"), _nowire("!(hh._stereo||hh._tdiv)"),
    _thit("hhmc._rel==0"),_bkghit("hhmc._rel!=0"),
    _crcan(0), _cpcan(0), _rcan(0), _fzcan(0), _corrcan(0), _hpcan1(0), _hpcan2(0), _hdcan(0), _t0can(0),
    _rsum(0), _hmva(0), _hres(0), _circ(0)
    { }

    void CenterRes();
    void CenterPos();
    void Radius();
    void RadiusSum();
    void PhiZ();
    void Correlations();
    void HitPos();
    void HitDist();
    void HitMVA();
    void HitRes();
    void Circle();
    void time();
    void save(const char* suffix=".png");

    TTree* _hdiag;
    TCut _helixOK;
    TCut _mchelixOK;
    TCut _hhout;
    TCut _stereo;
    TCut _tdiv;
    TCut _tdivonly;
    TCut _resphi;
    TCut _nowire;
    TCut _thit;
    TCut _bkghit;

    TCanvas *_crcan, *_cpcan, *_rcan, *_fzcan, *_corrcan, *_hpcan1, *_hpcan2, *_hdcan, *_t0can, 
	    *_rsum, *_hmva, *_hres, *_circ;
};

void HelixDiag::CenterRes() {
  TH2F* crcomp = new TH2F("crcomp","True Center Radius vs Reco;Reco Center Radius(mm);MC Center Radius(mm)",50,150.0,450.0,50,150.0,450.0);
  TH1F* crres = new TH1F("crres","Center Radius Resolution;#Delta R (mm)",100,-150.0,150.0);
  TH1F* ccrres = new TH1F("ccrres","Center Radius Resolution;#Delta R (mm)",100,-150.0,150.0);
  TH2F* cfcomp = new TH2F("cfcomp","Reco vs true Center #Phi;MC Center #phi;Reco Center #phi",50,-3.15,3.15,50,-3.15,3.15);
  TH1F* cfres = new TH1F("cfres","Center Radius #Phi Resolution;#Delta #phi",100,-0.4,0.4);
  crres->SetLineColor(kBlue);
  ccrres->SetLineColor(kRed);

  _hdiag->Project("crcomp","mch._rcent:rhel._rcent",_helixOK&&_mchelixOK);
  _hdiag->Project("crres","rhel._rcent-mch._rcent",_helixOK&&_mchelixOK);
  _hdiag->Project("cfcomp","rhel._fcent:mch._fcent",_helixOK&&_mchelixOK);
  _hdiag->Project("cfres","rhel._fcent-mch._fcent",_helixOK&&_mchelixOK);

  crcomp->FitSlicesY(0,0,-1,10);
  TH1D *crcomp_1 = (TH1D*)gDirectory->Get("crcomp_1");
  TH1D *crcomp_2 = (TH1D*)gDirectory->Get("crcomp_2");
  TH1D *crcomp_0 = (TH1D*)gDirectory->Get("crcomp_0");
  std::vector<Double_t> rreco, rmc, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(crcomp_0->GetBinContent(ibin+1)>0){
      rreco.push_back(crcomp_1->GetBinCenter(ibin+1));
      rmc.push_back(crcomp_1->GetBinContent(ibin+1));
      rerr.push_back(crcomp_2->GetBinContent(ibin+1)/sqrt(crcomp_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* me(0);
  if(rmc.size() > 2){
    me = new TGraphErrors(rmc.size(),rreco.data(),rmc.data(),rerr.data(),0);
    me->SetMarkerStyle(20);
    me->SetMarkerColor(kBlue);
    me->SetLineColor(kBlue);
    me->SetMarkerSize(0.5);
    //  rf->SetParameters(260.0,1.2,0.8);
    TF1* crf = new TF1("crf","[0]+[1]*x");
    crf->SetParameters(-100,1.5);
    TFitResultPtr crfit = me->Fit(crf,"S","same");
    char rcorr[100];
    //  snprintf(rcorr,100,"%f+(rhel._rcent>%f?%f:%f)*(rhel._rcent-%f)-mch._rcent",crfit->Parameter(0),crfit->Parameter(0),crfit->Parameter(1),crfit->Parameter(2),crfit->Parameter(0));
    snprintf(rcorr,100,"%f+rhel._rcent*%f-mch._rcent",crfit->Parameter(0),crfit->Parameter(1));
    cout << "projection = " << rcorr << endl;
    _hdiag->Project("ccrres",rcorr,_helixOK&&_mchelixOK);
  }
  _crcan = new TCanvas("crcan","Center Resolution",700,700);
  _crcan->Divide(2,2);
  _crcan->cd(1);
  crcomp->Draw();
  if(me != 0)me->Draw("LP");
//  TF1* rf = new TF1("rf","[0]+(x>[0]?[1]:[2])*(x-[0])");
  _crcan->cd(2);
  crres->Draw();
  ccrres->Draw("same");
  _crcan->cd(3);
  cfcomp->Draw();
  _crcan->cd(4);
  cfres->Draw();
}

void HelixDiag::CenterPos() {
  TH2F* rcpos = new TH2F("rcpos","Reco Center Position; x(mm); y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* mccpos = new TH2F("mccpos","MC Center Position; x(mm); y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* xcomp = new TH2F("xcomp","Reco vs true Center x;MC Center x (mm);Reco Center x (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* ycomp = new TH2F("ycomp","Reco vs true Center y;MC Center y (mm);Reco Center y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  _hdiag->Project("rcpos","rhel._rcent*sin(rhel._fcent):rhel._rcent*cos(rhel._fcent)",_helixOK&&_mchelixOK);
  _hdiag->Project("mccpos","mch._rcent*sin(mch._fcent):mch._rcent*cos(mch._fcent)",_helixOK&&_mchelixOK);
  _hdiag->Project("xcomp","rhel._rcent*cos(rhel._fcent):mch._rcent*cos(mch._fcent)",_helixOK&&_mchelixOK);
  _hdiag->Project("ycomp","rhel._rcent*sin(rhel._fcent):mch._rcent*sin(mch._fcent)",_helixOK&&_mchelixOK);

  _cpcan = new TCanvas("cpcan","Center Position",800,800);
  _cpcan->Divide(2,2);
  _cpcan->cd(1);
  rcpos->Draw();
  _cpcan->cd(2);
  mccpos->Draw();
  _cpcan->cd(3);
  xcomp->Draw();
  _cpcan->cd(4);
  ycomp->Draw();
}

void HelixDiag::Radius() {
  TH2F* rcomp = new TH2F("rcomp","True vs Reco Radius;Reco radius (mm); MC radius (mm)",50,200.0,350.0,50,200.0,350.0);
  TH1F* rres = new TH1F("rres","Radius resolution;reco - MC radius (mm)",100,-100,100);
  TH1F* rresc = new TH1F("rresc","Radius resolution;reco - MC radius (mm)",100,-100,100);
  rres->SetLineColor(kBlue);
  rresc->SetLineColor(kRed);
  _hdiag->Project("rcomp","mch._radius:rhel._radius",_helixOK&&_mchelixOK);
  _hdiag->Project("rres","rhel._radius-mch._radius",_helixOK&&_mchelixOK);
  
  rcomp->FitSlicesY(0,0,-1,20);
  TH1D *rcomp_1 = (TH1D*)gDirectory->Get("rcomp_1");
  TH1D *rcomp_2 = (TH1D*)gDirectory->Get("rcomp_2");
  TH1D *rcomp_0 = (TH1D*)gDirectory->Get("rcomp_0");
  std::vector<Double_t> rreco, rmc, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(rcomp_0->GetBinContent(ibin+1)>0){
      rreco.push_back(rcomp_1->GetBinCenter(ibin+1));
      rmc.push_back(rcomp_1->GetBinContent(ibin+1));
      rerr.push_back(rcomp_2->GetBinContent(ibin+1)/sqrt(rcomp_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* re(0);
  if(rmc.size() > 2){
    re = new TGraphErrors(rmc.size(),rreco.data(),rmc.data(),0,rerr.data());
    re->SetMarkerStyle(20);
    re->SetMarkerColor(kBlue);
    re->SetMarkerSize(0.5);
    re->SetLineColor(kBlue);
    //  rf->SetParameters(0.0,5.0,-0.0001);
    TF1* rf = new TF1("rf","[1]*(x-[0]) + pow(1.0+[1]*[1],1.5)*[2]*(x-[0])*(x-[0])+[3]");
    rf->SetParameters(260.0,1.0,-0.0001,260.0);
    TFitResultPtr rfit = re->Fit(rf,"S","same");
    char rcorr[120];
    //  snprintf(rcorr,100,"%f+(rhel._radius>%f?%f:%f)*(rhel._radius-%f)-mch._radius",rfit->Parameter(3),rfit->Parameter(0),rfit->Parameter(1),rfit->Parameter(2),rfit->Parameter(0));
    //  snprintf(rcorr,100,"%f+(rhel._radius-%f)*((rhel._radius>%f)?%f:%f)-mch._radius",rfit->Parameter(0),rfit->Parameter(0),rfit->Parameter(0),rfit->Parameter(1),rfit->Parameter(2));
    snprintf(rcorr,120,"%f+%f*(rhel._radius-%f)+pow(1.0+%f^2,1.5)*%f*(rhel._radius-%f)^2-mch._radius",rfit->Parameter(3),rfit->Parameter(1),rfit->Parameter(0),rfit->Parameter(1),rfit->Parameter(2),rfit->Parameter(0));
    cout << "projection = " << rcorr << endl;
    _hdiag->Project("rresc",rcorr,_helixOK&&_mchelixOK);
  }
  _rcan = new TCanvas("rcan","Radius",1000,500);
  _rcan->Divide(2,1);
  _rcan->cd(1);
  rcomp->Draw();
  if(re != 0)re->Draw("LP");
  //  TF1* rf = new TF1("rf","[3]+(x>[0]?[1]:[2])*(x-[0])");
  //  rf->SetParameters(240.0,0.4,0.8,290.0);
  //  TF1* rf = new TF1("rf","[0]+[1]*x+[2]*x*x");
  _rcan->cd(2);
  rres->Draw();
  rresc->Draw("same");
  TLegend* rrleg = new TLegend(0.6,0.7,0.9,0.9);
  rrleg->AddEntry(rres,"Uncorrectoed","L");
  rrleg->AddEntry(rresc,"Corrected","L");
  rrleg->Draw();
}

void HelixDiag::PhiZ() {
  TH2F* lcomp = new TH2F("lcomp","Reco vs true Lambda;MC #Lambda (mm);Reco #Lambda (mm)",50,125.0,325.0,50,125.0,325.0);
  TH1F* lres = new TH1F("lres","Lambda resolution;reco - MC lambda (mm)",100,-50,50);
  TH2F* fcomp = new TH2F("fcomp","Reco vs true #phiz0;MC #phiz0 (rad);Reco #fz0 (rad)",50,-3.15,3.15,50,-3.15,3.15);
  TH1F* fres = new TH1F("fres","#phiz0 resolution;Reco - MC fz0 (rad)",100,-0.4,0.4);

  _hdiag->Project("lcomp","rhel._lambda:mch._lambda",_helixOK&&_mchelixOK);
  _hdiag->Project("lres","rhel._lambda-mch._lambda",_helixOK&&_mchelixOK);
  _hdiag->Project("fcomp","rhel._fz0:mch._fz0",_helixOK&&_mchelixOK);
  _hdiag->Project("fres","rhel._fz0-mch._fz0",_helixOK&&_mchelixOK);

  _fzcan = new TCanvas("fzcan","Center Position",800,800);
  _fzcan->Divide(2,2);
  _fzcan->cd(1);
  lcomp->Draw();
  _fzcan->cd(2);
  lres->Draw();
  _fzcan->cd(3);
  fcomp->Draw();
  _fzcan->cd(4);
  fres->Draw();
}

void HelixDiag::Correlations() {
  TH2F* lrcorr = new TH2F("lrcorr","Reco - MC Lambda vs Radius;#Delta R (mm);#Delta #Lambda (mm)",50,-100.0,100.0,50,-100.0,100.0);
  TH2F* crrcorr = new TH2F("crrcorr","Reco - MC Center Radius vs Radius;#Delta R (mm);#Delta RC (mm)",50,-100.0,100.0,50,-100.0,100.0);
  TH2F* crlcorr = new TH2F("crlcorr","Reco - MC Center Radius vs Lambda;#Delta #Lambda (mm);#Delta RC (mm)",50,-100.0,100.0,50,-100.0,100.0);
  TH2F* cffcorr = new TH2F("cffcorr","Reco vs true center #phi vs #phiz0;#Delta #phi_{z0};#Delta #phi_{C}",50,-0.4,0.4,50,-0.4,0.4);
  _hdiag->Project("lrcorr","rhel._lambda-mch._lambda:rhel._radius-mch._radius",_helixOK&&_mchelixOK);
  _hdiag->Project("crrcorr","rhel._rcent-mch._rcent:rhel._radius-mch._radius",_helixOK&&_mchelixOK);
  _hdiag->Project("crlcorr","rhel._rcent-mch._rcent:rhel._lambda-mch._lambda",_helixOK&&_mchelixOK);
  _hdiag->Project("cffcorr","rhel._fcent-mch._fcent:rhel._fz0-mch._fz0",_helixOK&&_mchelixOK);

  crrcorr->FitSlicesY(0,0,-1,10);
  TH1D *crrcorr_1 = (TH1D*)gDirectory->Get("crrcorr_1");
  TH1D *crrcorr_2 = (TH1D*)gDirectory->Get("crrcorr_2");
  TH1D *crrcorr_0 = (TH1D*)gDirectory->Get("crrcorr_0");
  std::vector<Double_t> rcent, radius, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(crrcorr_0->GetBinContent(ibin+1)>0){
      radius.push_back(crrcorr_1->GetBinCenter(ibin+1));
      rcent.push_back(crrcorr_1->GetBinContent(ibin+1));
      rerr.push_back(crrcorr_2->GetBinContent(ibin+1)/sqrt(crrcorr_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* rrc = new TGraphErrors(radius.size(),radius.data(),rcent.data(),0,rerr.data());
  rrc->SetMarkerStyle(22);
  rrc->SetMarkerColor(kGreen);

 _corrcan = new TCanvas("corrcan","Corrrclations",800,800);
  _corrcan->Divide(2,2);
  _corrcan->cd(1);
  lrcorr->Draw();
  _corrcan->cd(2);
  crrcorr->Draw();
  rrc->Draw("same");
  rrc->Fit("pol1","","same");
  _corrcan->cd(3);
  crlcorr->Draw();
  _corrcan->cd(4);
  cffcorr->Draw();

}

void HelixDiag::Circle() {
  TH1F* hdrs = new TH1F("hdrs","Hit #Delta R",120,-150.0,150.0);
  TH1F* hdrb = new TH1F("hdrb","Hit #Delta R",120,-150.0,150.0);
  TH1F* hdrr = new TH1F("hdrr","Hit #Delta R",120,-150.0,150.0);
  hdrs->SetLineColor(kGreen);
  hdrb->SetLineColor(kRed);
  hdrr->SetLineColor(kBlue);
  hdrs->SetStats(0);
  hdrb->SetStats(0);
  hdrr->SetStats(0);

  TH1F* hrpulls = new TH1F("hrpulls","Hit #Delta R pull",120,-15.0,15.0);
  TH1F* hrpullb = new TH1F("hrpullb","Hit #Delta R pull",120,-15.0,15.0);
  TH1F* hrpullr = new TH1F("hrpullr","Hit #Delta R pull",120,-15.0,15.0);
  hrpulls->SetLineColor(kGreen);
  hrpullb->SetLineColor(kRed);
  hrpullr->SetLineColor(kBlue);
  hrpulls->SetStats(0);
  hrpullb->SetStats(0);
  hrpullr->SetStats(0);

  TH1F* hrwdots = new TH1F("hrwdots","Hit RW dot",200,-1.1,1.1);
  TH1F* hrwdotb = new TH1F("hrwdotb","Hit RW dot",200,-1.1,1.1);
  TH1F* hrwdotr = new TH1F("hrwdotr","Hit RW dot",200,-1.1,1.1);
  hrwdots->SetLineColor(kGreen);
  hrwdotb->SetLineColor(kRed);
  hrwdotr->SetLineColor(kBlue);
  hrwdots->SetStats(0);
  hrwdotb->SetStats(0);
  hrwdotr->SetStats(0);

  _hdiag->Project("hdrs","hh._hrho-rhel._radius",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hdrb","hh._hrho-rhel._radius",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hdrr","hh._hrho-rhel._radius",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hdiag->Project("hrpulls","(hh._hrho-rhel._radius)/max(hh._werr*hh._whdot,20.0)",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hrpullb","(hh._hrho-rhel._radius)/max(hh._werr*hh._whdot,20.0)",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hrpullr","(hh._hrho-rhel._radius)/max(hh._werr*hh._whdot,20.0)",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hdiag->Project("hrwdots","hh._whdot",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hrwdotb","hh._whdot",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hrwdotr","hh._whdot",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _circ = new TCanvas("circ","Circle Hit Resolution",800,800);
  _circ->Divide(2,2);
  _circ->cd(1);
  gPad->SetLogy();
  hdrs->Draw();
  hdrb->Draw("same");
  hdrr->Draw("same");
  _circ->cd(2);
  gPad->SetLogy();
  hrpulls->Draw();
  hrpullb->Draw("same");
  hrpullr->Draw("same");
  _circ->cd(3);
  hrwdots->Draw();
  hrwdotb->Draw("same");
  hrwdotr->Draw("same");

  TLegend* lcirc = new TLegend(0.2,0.7,0.5,0.9);
  lcirc->AddEntry(hdrs,"Direct","L");
  lcirc->AddEntry(hdrb,"Background","L");
  lcirc->AddEntry(hdrr,"Related","L");
  lcirc->Draw();

}

void HelixDiag::HitRes() {
  TH1F* hdphis = new TH1F("hdphis","Hit #Delta phi",120,-2,2);
  TH1F* hdphib = new TH1F("hdphib","Hit #Delta phi",120,-2,2);
  TH1F* hdphir = new TH1F("hdphir","Hit #Delta phi",120,-2,2);
  hdphis->SetLineColor(kGreen);
  hdphib->SetLineColor(kRed);
  hdphir->SetLineColor(kBlue);
  hdphis->SetStats(0);
  hdphib->SetStats(0);
  hdphir->SetStats(0);

  TH1F* hdwires = new TH1F("hdwires","Hit #Delta wire;#Delta dist (mm)",120,-250.0,250.0);
  TH1F* hdwireb = new TH1F("hdwireb","Hit #Delta wire;#Delta dist (mm)",120,-250.0,250.0);
  TH1F* hdwirer = new TH1F("hdwirer","Hit #Delta wire;#Delta dist (mm)",120,-250.0,250.0);
  hdwires->SetLineColor(kGreen);
  hdwireb->SetLineColor(kRed);
  hdwirer->SetLineColor(kBlue);
  hdwires->SetStats(0);
  hdwireb->SetStats(0);
  hdwirer->SetStats(0);

  TH1F* hdtranss = new TH1F("hdtranss","Hit #Delta trans;#Delta dist (mm)",120,-100.0,100.0);
  TH1F* hdtransb = new TH1F("hdtransb","Hit #Delta trans;#Delta dist (mm)",120,-100.0,100.0);
  TH1F* hdtransr = new TH1F("hdtransr","Hit #Delta trans;#Delta dist (mm)",120,-100.0,100.0);
  hdtranss->SetLineColor(kGreen);
  hdtransb->SetLineColor(kRed);
  hdtransr->SetLineColor(kBlue);
  hdtranss->SetStats(0);
  hdtransb->SetStats(0);
  hdtransr->SetStats(0);

  TH1F* hchisqs = new TH1F("hchisqs","Hit chisq",300,0.0,150.0);
  TH1F* hchisqb = new TH1F("hchisqb","Hit chisq",300,0.0,150.0);
  TH1F* hchisqr = new TH1F("hchisqr","Hit chisq",300,0.0,150.0);
  hchisqs->SetLineColor(kGreen);
  hchisqb->SetLineColor(kRed);
  hchisqr->SetLineColor(kBlue);
  hchisqs->SetStats(0);
  hchisqb->SetStats(0);
  hchisqr->SetStats(0);

  _hdiag->Project("hdphis","atan2(hh._hhpos.dy,hh._hhpos.dx)-rhel._fcent",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hdphib","atan2(hh._hhpos.dy,hh._hhpos.dx)-rhel._fcent",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hdphir","atan2(hh._hhpos.dy,hh._hhpos.dx)-rhel._fcent",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hdiag->Project("hdwires","hh._dwire",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hdwireb","hh._dwire",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hdwirer","hh._dwire",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hdiag->Project("hdtranss","hh._dtrans",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hdtransb","hh._dtrans",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hdtransr","hh._dtrans",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hdiag->Project("hchisqs","(hh._dwire/hh._wres)^2+(hh._dtrans/hh._wtres)^2",_helixOK&&_mchelixOK&&"hhmc._rel==0");
  _hdiag->Project("hchisqb","(hh._dwire/hh._wres)^2+(hh._dtrans/hh._wtres)^2",_helixOK&&_mchelixOK&&"hhmc._rel<0");
  _hdiag->Project("hchisqr","(hh._dwire/hh._wres)^2+(hh._dtrans/hh._wtres)^2",_helixOK&&_mchelixOK&&"hhmc._rel>0");

  _hres = new TCanvas("hres","Hit Resolution",800,800);
  _hres->Divide(2,2);
  _hres->cd(1);
  hdphis->Draw();
  hdphib->Draw("same");
  hdphir->Draw("same");
  _hres->cd(2);
  gPad->SetLogy();
  hdwires->Draw();
  hdwireb->Draw("same");
  hdwirer->Draw("same");
  _hres->cd(3);
  gPad->SetLogy();
  hdtranss->Draw();
  hdtransb->Draw("same");
  hdtransr->Draw("same");
  _hres->cd(4);
  gPad->SetLogy();
  hchisqs->Draw();
  hchisqb->Draw("same");
  hchisqr->Draw("same");
  TLegend* lhres = new TLegend(0.2,0.7,0.5,0.9);
  lhres->AddEntry(hdphis,"Direct","L");
  lhres->AddEntry(hdphib,"Background","L");
  lhres->AddEntry(hdphir,"Related","L");
  lhres->Draw();
}

void HelixDiag::HitMVA() {
  TH1F* hmvanos = new TH1F("hmvanos","Hit MVA, Not Outlier",120,-0.1,1.1);
  TH1F* hmvanob = new TH1F("hmvanob","Hit MVA, Not Outlier",120,-0.1,1.1);
  TH1F* hmvanor = new TH1F("hmvanor","Hit MVA, Not Outlier",120,-0.1,1.1);
  hmvanos->SetLineColor(kGreen);
  hmvanob->SetLineColor(kRed);
  hmvanor->SetLineColor(kBlue);
  TH1F* hmvaots = new TH1F("hmvaots","Hit MVA, Outlier",120,-0.1,1.1);
  TH1F* hmvaotb = new TH1F("hmvaotb","Hit MVA, Outlier",120,-0.1,1.1);
  TH1F* hmvaotr = new TH1F("hmvaotr","Hit MVA, Outlier",120,-0.1,1.1);
  hmvaots->SetLineColor(kGreen);
  hmvaotb->SetLineColor(kRed);
  hmvaotr->SetLineColor(kBlue);

  hmvanos->SetStats(0);
  hmvanob->SetStats(0);
  hmvanor->SetStats(0);
  hmvaots->SetStats(0);
  hmvaotb->SetStats(0);
  hmvaotr->SetStats(0);

  _hdiag->Project("hmvanos","hh._hqual",_helixOK&&_mchelixOK&&(!_hhout)&&"hhmc._rel==0");
  _hdiag->Project("hmvanob","hh._hqual",_helixOK&&_mchelixOK&&(!_hhout)&&"hhmc._rel<0");
  _hdiag->Project("hmvanor","hh._hqual",_helixOK&&_mchelixOK&&(!_hhout)&&"hhmc._rel>0");

  _hdiag->Project("hmvaots","hh._hqual",_helixOK&&_mchelixOK&&_hhout&&"hhmc._rel==0");
  _hdiag->Project("hmvaotb","hh._hqual",_helixOK&&_mchelixOK&&_hhout&&"hhmc._rel<0");
  _hdiag->Project("hmvaotr","hh._hqual",_helixOK&&_mchelixOK&&_hhout&&"hhmc._rel>0");

  _hmva = new TCanvas("hmva","Hit MVA", 800,600);
  _hmva->Divide(2,1);
  _hmva->cd(1);
  hmvanos->Draw();
  hmvanob->Draw("same");
  hmvanor->Draw("same");
  _hmva->cd(2);
  hmvaotb->Draw();
  hmvaots->Draw("same");
  hmvaotr->Draw("same");

  TLegend* mval = new TLegend(0.2,0.7,0.5,0.9);
  mval->AddEntry(hmvanos,"Direct","L");
  mval->AddEntry(hmvanob,"Background","L");
  mval->AddEntry(hmvanor,"Related","L");
  mval->Draw();
}

void HelixDiag::HitPos() {
  TH2F* hrcomp = new TH2F("hrcomp","Reco vs true Hit radius",100,350.0,700.0,100,350.0,700.0);
  TH2F* hfcomp = new TH2F("hfcomp","Reco vs true Hit #phi",100,-3.15,3.15,100,-3.15,3.15);
  TH2F* hercomp = new TH2F("hercomp","Expected vs true Hit radius",100,350.0,700.0,100,350.0,700.0);
  TH2F* hefcomp = new TH2F("hefcomp","Expected vs true Hit #phi",100,-3.15,3.15,100,-3.15,3.15);
  _hdiag->Project("hrcomp","sqrt(hh._hhpos.dy^2+hh._hhpos.dx^2):sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("hfcomp","atan2(hh._hhpos.dy,hh._hhpos.dx):atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("hercomp","sqrt(hh._hpos.dy^2+hh._hpos.dx^2):sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("hefcomp","atan2(hh._hpos.dy,hh._hpos.dx):atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mchelixOK&&!_hhout);

 _hpcan1 = new TCanvas("hpcan1","Hit Positions",800,800);
  _hpcan1->Divide(2,2);
  _hpcan1->cd(1);
  hrcomp->Draw("colorz");
  _hpcan1->cd(2);
  hfcomp->Draw("colorz");
  _hpcan1->cd(3);
  hercomp->Draw("colorz");
  _hpcan1->cd(4);
  hefcomp->Draw("colorz");

  TH1F* hrres = new TH1F("hrres","Reco vs true Hit radius",100,-100.0,100.0);
  TH1F* hfres = new TH1F("hfres","Reco vs true Hit #phi",100,-0.4,0.4);
  TH1F* herres = new TH1F("herres","Expected vs true Hit radius",100,-100.0,100.0);
  TH1F* hefres = new TH1F("hefres","Expected vs true Hit #phi",100,-0.4,0.4);
  _hdiag->Project("hrres","sqrt(hh._hhpos.dy^2+hh._hhpos.dx^2)-sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("hfres","atan2(hh._hhpos.dy,hh._hhpos.dx)-atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("herres","sqrt(hh._hpos.dy^2+hh._hpos.dx^2)-sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mchelixOK&&!_hhout);
  _hdiag->Project("hefres","atan2(hh._hpos.dy,hh._hpos.dx)-atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mchelixOK&&!_hhout);

 _hpcan2 = new TCanvas("hpcan2","Hit Position Resolution",800,800);
  _hpcan2->Divide(2,2);
  _hpcan2->cd(1);
  hrres->Draw();
  _hpcan2->cd(2);
  hfres->Draw();
  _hpcan2->cd(3);
  herres->Draw();
  _hpcan2->cd(4);
  hefres->Draw();

}

void HelixDiag::HitDist() {
  TH1F* hdwires = new TH1F("hdwires","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  TH1F* hdwiret = new TH1F("hdwiret","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  TH1F* hdwireo = new TH1F("hdwireo","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  hdwires->SetLineColor(kBlue);
  hdwiret->SetLineColor(kGreen);
  hdwireo->SetLineColor(kBlack);
  hdwires->SetStats(0);
  hdwiret->SetStats(0);
  hdwireo->SetStats(0);
  _hdiag->Project("hdwires","hh._dwire",_helixOK&&_mchelixOK&&_thit&&(!_hhout)&&_stereo);
  _hdiag->Project("hdwiret","hh._dwire",_helixOK&&_mchelixOK&&_thit&&(!_hhout)&&_tdivonly);
  _hdiag->Project("hdwireo","hh._dwire",_helixOK&&_mchelixOK&&_thit&&_hhout);

  TH1F* hdtranss = new TH1F("hdtranss","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtranst = new TH1F("hdtranst","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtranso = new TH1F("hdtranso","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  hdtranss->SetLineColor(kBlue);
  hdtranst->SetLineColor(kGreen);
  hdtranso->SetLineColor(kBlack);
  hdtranss->SetStats(0);
  hdtranst->SetStats(0);
  hdtranso->SetStats(0);
  _hdiag->Project("hdtranss","hh._dtrans",_helixOK&&_mchelixOK&&_thit&&(!_hhout)&&_stereo);
  _hdiag->Project("hdtranst","hh._dtrans",_helixOK&&_mchelixOK&&_thit&&(!_hhout)&&_tdivonly);
  _hdiag->Project("hdtranso","hh._dtrans",_helixOK&&_mchelixOK&&_thit&&_hhout);

  TLegend* dleg = new TLegend(0.6,0.7,0.9,0.9);
  dleg->AddEntry(hdwires,"Stereo Hits","L");
  dleg->AddEntry(hdwiret,"TimeDiv Hits","L");
  dleg->AddEntry(hdwireo,"Outlier Hits","L");
  _hdcan = new TCanvas("hdcan","Hit Dists",800,800);
  _hdcan->Divide(2,2);
  _hdcan->cd(1);
  hdwires->Draw();
  hdwiret->Draw("same");
  hdwireo->Draw("same");
  dleg->Draw();
  _hdcan->cd(2);
  hdtranss->Draw();
  hdtranst->Draw("same");
  hdtranso->Draw("same");
}

void HelixDiag::time() {
  TH1F* tct0 = new TH1F("tct0","T0 Resolution;t0_{reco}-t0_{MC} (ns)",100,-20,20);
  TH1F* ht0 = new TH1F("ht0","T0 Resolution;t0_{reco}-t0_{MC} (ns)",100,-20,20);
  tct0->SetLineColor(kBlue);
  ht0->SetLineColor(kRed);
  _hdiag->Project("tct0","tct0-mct0%1695",_helixOK&&_mchelixOK);
  _hdiag->Project("ht0","ht0-mct0%1695",_helixOK&&_mchelixOK);

  TH1F* thdt = new TH1F("thdt","Hit #Delta t;t_{hit}-t0",100,-25,75);
  TH1F* bhdt = new TH1F("bhdt","Hit #Delta t;t_{hit}-t0",100,-25,75);
  _hdiag->Project("thdt","hh._dt",_helixOK&&_mchelixOK&&_thit);
  _hdiag->Project("bhdt","hh._dt",_helixOK&&_mchelixOK&&_bkghit);
  thdt->SetLineColor(kGreen);
  bhdt->SetLineColor(kCyan);

  _t0can = new TCanvas("t0can","T0",600,600);
  _t0can->Divide(2,2);
  _t0can->cd(1);
  ht0->Draw();
  tct0->Draw("same");
  TLegend* t0leg = new TLegend(0.1,0.6,0.3,0.9);
  t0leg->AddEntry(tct0,"Time Cluster","L");
  t0leg->AddEntry(ht0,"Helix","L");
  t0leg->Draw();
  _t0can->cd(2);
  thdt->Draw();
  bhdt->Draw("same");
}

void HelixDiag::RadiusSum() {
  TH2F* drcomp = new TH2F("drcomp","Reco vs True Radius Diff;MC #Delta (mm);Reco #Delta R (mm)",50,-150.0,150.0,50,-150.0,150.0);
  TH2F* srcomp = new TH2F("srcomp","Reco vs True Radius Sum;MC #Sigma (mm);Reco #Sigma R (mm)",50,400.0,700.0,50,400.0,700.0);

  TH1F* drres = new TH1F("drres","Radius Diff Resolution;#Delta R (mm)",100,-150.0,150.0);
  TH1F* srres = new TH1F("srres","Radius Sum Resolution;#Delta R (mm)",100,-150.0,150.0);

  _hdiag->Project("drcomp","rhel._radius-rhel._rcent:mch._radius-mch._rcent",_helixOK&&_mchelixOK);
  _hdiag->Project("drres","(rhel._radius-rhel._rcent)-(mch._radius-mch._rcent)",_helixOK&&_mchelixOK);

  _hdiag->Project("srcomp","rhel._radius+rhel._rcent:mch._radius+mch._rcent",_helixOK&&_mchelixOK);
  _hdiag->Project("srres","(rhel._radius+rhel._rcent)-(mch._radius+mch._rcent)",_helixOK&&_mchelixOK);

  srcomp->FitSlicesX(0,0,-1,10);
  TH1D *srcomp_1 = (TH1D*)gDirectory->Get("srcomp_1");
  TH1D *srcomp_2 = (TH1D*)gDirectory->Get("srcomp_2");
  TH1D *srcomp_0 = (TH1D*)gDirectory->Get("srcomp_0");
  std::vector<Double_t> rreco, rmc, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(srcomp_0->GetBinContent(ibin+1)>0){
      rmc.push_back(srcomp_1->GetBinCenter(ibin+1));
      rreco.push_back(srcomp_1->GetBinContent(ibin+1));
      rerr.push_back(srcomp_2->GetBinContent(ibin+1)/sqrt(srcomp_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* me = new TGraphErrors(rmc.size(),rmc.data(),rreco.data(),0,rerr.data());
  me->SetMarkerStyle(20);
  me->SetMarkerColor(kBlue);
  me->SetLineColor(kBlue);
  me->SetMarkerSize(0.5);
  TF1* crf = new TF1("crf","[0]+[1]*x");
  crf->SetParameters(-100,1.5);
  TFitResultPtr crfit = me->Fit(crf,"S","same");

  _rsum = new TCanvas("rsum","rsum",700,700);
  _rsum->Divide(2,2);
  _rsum->cd(1);
  drcomp->Draw();
  _rsum->cd(2);
  drres->Draw();
  _rsum->cd(3);
  srcomp->Draw();
  me->Draw("PL");

  _rsum->cd(4);
  srres->Draw();
}

void HelixDiag::save(const char* suffix) {
  string ss(suffix);
  string cfname;
  if(_crcan){
    cfname = string(_crcan->GetName())+ss;
    _crcan->SaveAs(cfname.c_str());
  }
  if(_cpcan){
    cfname = string(_cpcan->GetName())+ss;
    _cpcan->SaveAs(cfname.c_str());
  }
  if(_rcan){
    cfname = string(_rcan->GetName())+ss;
    _rcan->SaveAs(cfname.c_str());
  }
  if(_fzcan){
    cfname = string(_fzcan->GetName())+ss;
    _fzcan->SaveAs(cfname.c_str());
  }
  if(_corrcan){
    cfname = string(_corrcan->GetName())+ss;
    _corrcan->SaveAs(cfname.c_str());
  }
  if(_hpcan1){
    cfname = string(_hpcan1->GetName())+ss;
    _hpcan1->SaveAs(cfname.c_str());
  }
  if(_hpcan2){
    cfname = string(_hpcan2->GetName())+ss;
    _hpcan2->SaveAs(cfname.c_str());
  }
  if(_hdcan){
    cfname = string(_hdcan->GetName())+ss;
    _hdcan->SaveAs(cfname.c_str());
  }
}

