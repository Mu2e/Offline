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
#include "TGraphErrors.h"

class HelixDiag  {
  public:
    HelixDiag(TTree* hdiag,TCut* mcsel=0) : _hdiag(hdiag), _helixOK("helixOK"), _mchelixOK("mchelixOK"),
    _outlier("hh._outlier"), _stereo("hh._stereo"), _tdiv("hh._tdiv"),
    _tdivonly("hh._tdiv && !hh._stereo"), _resphi("hh._resphi"), _nowire("!(hh._stereo||hh._tdiv)"),
    _crcan(0), _cpcan(0), _rcan(0), _fzcan(0), _corrcan(0), _hpcan1(0), _hpcan2(0), _hdcan(0)
    { if(mcsel != 0)_mcsel = *mcsel;}

    void CenterRes();
    void CenterPos();
    void Radius();
    void PhiZ();
    void Correlations();
    void HitPos();
    void HitDist();
    void save(const char* suffix=".png");

    TTree* _hdiag;
    TCut _mcsel;
    TCut _helixOK;
    TCut _mchelixOK;
    TCut _outlier;
    TCut _stereo;
    TCut _tdiv;
    TCut _tdivonly;
    TCut _resphi;
    TCut _nowire;

    TCanvas *_crcan, *_cpcan, *_rcan, *_fzcan, *_corrcan, *_hpcan1, *_hpcan2, *_hdcan;
};

void HelixDiag::CenterRes() {
  TH2F* crcomp = new TH2F("crcomp","Reco vs true Center Radius;MC Center Radius(mm);Reco Center Radius(mm)",50,150.0,450.0,50,150.0,450.0);
  TH1F* crres = new TH1F("crres","Center Radius Resolution;#Delta R (mm)",100,-150.0,150.0);
  TH2F* cfcomp = new TH2F("cfcomp","Reco vs true Center #Phi;MC Center #phi;Reco Center #phi",50,-3.15,3.15,50,-3.15,3.15);
  TH1F* cfres = new TH1F("cfres","Center Radius #Phi Resolution;#Delta #phi",100,-0.4,0.4);

  _hdiag->Project("crcomp","rhel._rcent:mch._rcent",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("crres","rhel._rcent-mch._rcent",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("cfcomp","rhel._fcent:mch._fcent",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("cfres","rhel._fcent-mch._fcent",_helixOK&&_mcsel&&_mchelixOK);

  crcomp->FitSlicesY(0,0,-1,10);
  TH1D *crcomp_1 = (TH1D*)gDirectory->Get("crcomp_1");
  TH1D *crcomp_2 = (TH1D*)gDirectory->Get("crcomp_2");
  TH1D *crcomp_0 = (TH1D*)gDirectory->Get("crcomp_0");
  std::vector<Double_t> rreco, rmc, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(crcomp_0->GetBinContent(ibin+1)>0){
      rmc.push_back(crcomp_1->GetBinCenter(ibin+1));
      rreco.push_back(crcomp_1->GetBinContent(ibin+1));
      rerr.push_back(crcomp_2->GetBinContent(ibin+1)/sqrt(crcomp_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* me = new TGraphErrors(rmc.size(),rmc.data(),rreco.data(),0,rerr.data());
  me->SetMarkerStyle(22);
  me->SetMarkerColor(kGreen);

  _crcan = new TCanvas("crcan","Center Resolution",800,800);
  _crcan->Divide(2,2);
  _crcan->cd(1);
  crcomp->Draw();
  me->Draw("same");
  me->Fit("pol1","","same");
  _crcan->cd(2);
  crres->Draw();
  _crcan->cd(3);
  cfcomp->Draw();
  _crcan->cd(4);
  cfres->Draw();
}

void HelixDiag::CenterPos() {
  TH2F* rcpos = new TH2F("rcpos","Reco Center Position; x(mm) y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* mccpos = new TH2F("mccpos","MC Center Position; x(mm) y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* xcomp = new TH2F("xcomp","Reco vs true Center x;Reco Center x (mm);MC Center x (mm)",50,-450.0,450.0,50,-450.0,450.0);
  TH2F* ycomp = new TH2F("ycomp","Reco vs true Center y;Reco Center y (mm);MC Center y (mm)",50,-450.0,450.0,50,-450.0,450.0);
  _hdiag->Project("rcpos","rhel._rcent*sin(rhel._fcent):rhel._rcent*cos(rhel._fcent)",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("mccpos","mch._rcent*sin(mch._fcent):mch._rcent*cos(mch._fcent)",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("xcomp","rhel._rcent*cos(rhel._fcent):mch._rcent*cos(mch._fcent)",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("ycomp","rhel._rcent*sin(rhel._fcent):mch._rcent*sin(mch._fcent)",_helixOK&&_mcsel&&_mchelixOK);

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
  TH2F* rcomp = new TH2F("rcomp","Reco vs true Radius;MC radius (mm); Reco radius (mm)",50,200.0,350.0,50,200.0,350.0);
  TH1F* rres = new TH1F("rres","Radius resolution;reco - MC radius (mm)",100,-100,100);
  _hdiag->Project("rcomp","rhel._radius:mch._radius",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("rres","rhel._radius-mch._radius",_helixOK&&_mcsel&&_mchelixOK);
  
  rcomp->FitSlicesY(0,0,-1,10);
  TH1D *rcomp_1 = (TH1D*)gDirectory->Get("rcomp_1");
  TH1D *rcomp_2 = (TH1D*)gDirectory->Get("rcomp_2");
  TH1D *rcomp_0 = (TH1D*)gDirectory->Get("rcomp_0");
  std::vector<Double_t> rreco, rmc, rerr;
  for(size_t ibin=0;ibin<50;++ibin){
    if(rcomp_0->GetBinContent(ibin+1)>0){
      rmc.push_back(rcomp_1->GetBinCenter(ibin+1));
      rreco.push_back(rcomp_1->GetBinContent(ibin+1));
      rerr.push_back(rcomp_2->GetBinContent(ibin+1)/sqrt(rcomp_0->GetBinContent(ibin+1)));
    }
  }
  TGraphErrors* re = new TGraphErrors(rmc.size(),rmc.data(),rreco.data(),0,rerr.data());
  re->SetMarkerStyle(22);
  re->SetMarkerColor(kGreen);

  _rcan = new TCanvas("rcan","Radius",800,600);
  _rcan->Divide(2,1);
  _rcan->cd(1);
  rcomp->Draw();
  re->Draw("same");
  re->Fit("pol1","","same");
  _rcan->cd(2);
  rres->Draw();


}

void HelixDiag::PhiZ() {
  TH2F* lcomp = new TH2F("lcomp","Reco vs true Lambda;MC #Lambda (mm);Reco #Lambda (mm)",50,125.0,325.0,50,125.0,325.0);
  TH1F* lres = new TH1F("lres","Lambda resolution;reco - MC lambda (mm)",100,-50,50);
  TH2F* fcomp = new TH2F("fcomp","Reco vs true #phiz0;MC #phiz0 (rad);Reco #fz0 (rad)",50,-3.15,3.15,50,-3.15,3.15);
  TH1F* fres = new TH1F("fres","#phiz0 resolution;Reco - MC fz0 (rad)",100,-0.4,0.4);

  _hdiag->Project("lcomp","rhel._lambda:mch._lambda",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("lres","rhel._lambda-mch._lambda",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("fcomp","rhel._fz0:mch._fz0",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("fres","rhel._fz0-mch._fz0",_helixOK&&_mcsel&&_mchelixOK);

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
  _hdiag->Project("lrcorr","rhel._lambda-mch._lambda:rhel._radius-mch._radius",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("crrcorr","rhel._rcent-mch._rcent:rhel._radius-mch._radius",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("crlcorr","rhel._rcent-mch._rcent:rhel._lambda-mch._lambda",_helixOK&&_mcsel&&_mchelixOK);
  _hdiag->Project("cffcorr","rhel._fcent-mch._fcent:rhel._fz0-mch._fz0",_helixOK&&_mcsel&&_mchelixOK);

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

void HelixDiag::HitPos() {
  TH2F* hrcomp = new TH2F("hrcomp","Reco vs true Hit radius",100,350.0,700.0,100,350.0,700.0);
  TH2F* hfcomp = new TH2F("hfcomp","Reco vs true Hit #phi",100,-3.15,3.15,100,-3.15,3.15);
  TH2F* hercomp = new TH2F("hercomp","Expected vs true Hit radius",100,350.0,700.0,100,350.0,700.0);
  TH2F* hefcomp = new TH2F("hefcomp","Expected vs true Hit #phi",100,-3.15,3.15,100,-3.15,3.15);
  _hdiag->Project("hrcomp","sqrt(hh._hhpos.dy^2+hh._hhpos.dx^2):sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("hfcomp","atan2(hh._hhpos.dy,hh._hhpos.dx):atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("hercomp","sqrt(hh._hpos.dy^2+hh._hpos.dx^2):sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("hefcomp","atan2(hh._hpos.dy,hh._hpos.dx):atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);

 _hpcan1 = new TCanvas("hpcan1","Hiti Positions",800,800);
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
  _hdiag->Project("hrres","sqrt(hh._hhpos.dy^2+hh._hhpos.dx^2)-sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("hfres","atan2(hh._hhpos.dy,hh._hhpos.dx)-atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("herres","sqrt(hh._hpos.dy^2+hh._hpos.dx^2)-sqrt(hhmc._hpos.dy^2+hhmc._hpos.dx^2)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);
  _hdiag->Project("hefres","atan2(hh._hpos.dy,hh._hpos.dx)-atan2(hhmc._hpos.dy,hhmc._hpos.dx)",_helixOK&&_mcsel&&_mchelixOK&&!_outlier);

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
  TH1F* hdwiren = new TH1F("hdwiren","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  TH1F* hdwireo = new TH1F("hdwireo","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  TH1F* hdwiref = new TH1F("hdwiref","Hit Wire Dist to Helix Center;D_{wire} (mm)",100,-500.0,500.0);
  hdwires->SetLineColor(kBlue);
  hdwiret->SetLineColor(kGreen);
  hdwiren->SetLineColor(kRed);
  hdwireo->SetLineColor(kBlack);
  hdwiref->SetLineColor(kCyan);
  hdwires->SetStats(0);
  hdwiret->SetStats(0);
  hdwiren->SetStats(0);
  hdwireo->SetStats(0);
  hdwiref->SetStats(0);
  _hdiag->Project("hdwires","hh._dwire",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&_stereo);
  _hdiag->Project("hdwiret","hh._dwire",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&_tdivonly);
  _hdiag->Project("hdwiren","hh._dwire",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&(!_tdiv)&&(!_stereo));
  _hdiag->Project("hdwireo","hh._dwire",_helixOK&&_mcsel&&_mchelixOK&&_outlier);
  _hdiag->Project("hdwiref","hh._dwire",!_helixOK&&_mcsel&&_mchelixOK&&!_outlier);

  TH1F* hdtranss = new TH1F("hdtranss","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtranst = new TH1F("hdtranst","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtransn = new TH1F("hdtransn","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtranso = new TH1F("hdtranso","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  TH1F* hdtransf = new TH1F("hdtransf","Hit Trans Dist to Helix Center;D_{trans} (mm)",100,-250.0,250.0);
  hdtranss->SetLineColor(kBlue);
  hdtranst->SetLineColor(kGreen);
  hdtransn->SetLineColor(kRed);
  hdtranso->SetLineColor(kBlack);
  hdtransf->SetLineColor(kCyan);
  hdtranss->SetStats(0);
  hdtranst->SetStats(0);
  hdtransn->SetStats(0);
  hdtranso->SetStats(0);
  hdtransf->SetStats(0);
  _hdiag->Project("hdtranss","hh._dtrans",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&_stereo);
  _hdiag->Project("hdtranst","hh._dtrans",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&_tdivonly);
  _hdiag->Project("hdtransn","hh._dtrans",_helixOK&&_mcsel&&_mchelixOK&&(!_outlier)&&(!_tdiv)&&(!_stereo));
  _hdiag->Project("hdtranso","hh._dtrans",_helixOK&&_mcsel&&_mchelixOK&&_outlier);
  _hdiag->Project("hdtransf","hh._dtrans",!_helixOK&&_mcsel&&_mchelixOK&&!_outlier);

  TLegend* dleg = new TLegend(0.6,0.7,0.9,0.9);
  dleg->AddEntry(hdwires,"Stereo Hits","L");
  dleg->AddEntry(hdwiret,"TimeDiv Hits","L");
  dleg->AddEntry(hdwiren,"No WireInfo Hits","L");
  dleg->AddEntry(hdwireo,"Outlier Hits","L");
  dleg->AddEntry(hdwiref,"Failed Fit","L");
  _hdcan = new TCanvas("hdcan","Hit Dists",800,800);
  _hdcan->Divide(2,2);
  _hdcan->cd(1);
  hdwires->Draw();
  hdwiret->Draw("same");
  hdwiren->Draw("same");
  hdwireo->Draw("same");
  hdwiref->Draw("same");
  dleg->Draw();
  _hdcan->cd(2);
  hdtranss->Draw();
  hdtranst->Draw("same");
  hdtransn->Draw("same");
  hdtranso->Draw("same");
  hdtransf->Draw("same");
 
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

