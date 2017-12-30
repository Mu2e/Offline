#define TDCalib_cxx
#include "TDCalib.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

Double_t ALine(Double_t *x, Double_t *par) {
  if(x[0] < par[0])
    return par[1];
  else
    return par[1]+(x[0]-par[0])*par[2];
}

Double_t DLine(Double_t *x, Double_t *par) {
  if(x[0] > par[0])
    return par[1]+(x[0]-par[0])*par[3];
  else
    return par[1]+(x[0]-par[0])*par[2];
}
 

void TDCalib::MakeHists(){
  _tdevec.clear();
  _tdpevec.clear();
  _evec.clear();
  _resvec.clear();
  char name[50];
  char profname[50];
  char ename[50];
  char rname[50];
  char title[80];
  char etitle[80];
  char rtitle[80];
  float emin = _emin;
  float emax = emin+_ebin;
  while(emax < _emax){
  // effective velocity (slope)
    snprintf(name,50,"td%3.1f",emin);
    snprintf(profname,50,"tdp%3.1f",emin);
    snprintf(title,80,"True wire position vs #Delta t, %3.1fKeV < edep < %3.1fKeV",emin,emax);
    TH2F* hist = new TH2F(name,title,50,-10.0,10.0,40,-600.0,600.0);
    _tdevec.push_back(hist);
    TProfile* prof = new TProfile(profname,title,50,-10.0,10.0,-600.0,600.0);
    _tdpevec.push_back(prof);
    // energy
    snprintf(ename,50,"edep%3.1f",emin);
    snprintf(etitle,80,"EDep, %3.1fKeV < edep < %3.1fKeV:EDep (KeV)",emin,emax);
    TH1F* edep = new TH1F(ename,etitle,100,0.0,7.0);
    _evec.push_back(edep);
    // resolution

    float rcore = 17.8 + 46.5*exp(-emin/1.7);
    snprintf(rname,50,"res%3.1f",emin);
    snprintf(rtitle,80,"TD Resolution vs wire position, %3.1fKeV < edep < %3.1fKeV;Wire Position(mm);Reco - mc wlen (mm)",emin,emax);
    TH2F* rhist = new TH2F(rname,rtitle,15,0.0,450.0,50,-6*rcore,6*rcore);
    _resvec.push_back(rhist);
      // next bin
    emin += _ebin;
    emax += _ebin;
   }
   _pulledep = new TH2F("pulledep","TD Pull vs Reco EDep;EDep (KeV);TD Pull",20,0.0,7.0,50,-10.0,10.0);
   _pullwlen = new TH2F("pullwlen","TD Pull vs Reco Wire Pos;Wire Pos (mm);TD Pull",20,0.0,650.0,50,-10.0,10.0);
   _pulltedep = new TH2F("pulltedep","TD Pull vs True EDep;MC EDep (KeV);TD Pull",20,0.0,7.0,50,-10.0,10.0);
   _pulltwlen = new TH2F("pulltwlen","TD Pull vs True Wire Pos;MC Wire Pos (mm);TD Pull",20,0.0,650.0,50,-10.0,10.0);
   _pullsfrac = new TH2F("pullsfrac","TD Pull vs Reco Wire Fraction;Wire pos/straw 1/2 length;TD Pull",20,0.0,1.01,50,-10.0,10.0);
   _pulltsfrac = new TH2F("pulltsfrac","TD Pull vs True Wire Fraction;MC Wire pos/straw 1/2 length;TD Pull",20,0.0,1.01,50,-10.0,10.0);
   _ceres = new TH1F("ceres","TD Ce Resolution;Res (mm)",100,-150,150);
   _pres = new TH1F("pres","TD Proton Resolution;Res (mm)",100,-150,150);
   _cepull = new TH1F("cepull","TD Ce Pull",100,-10,10);
   _ppull = new TH1F("ppull","TD Proton Pull",100,-10,10);
}

void TDCalib::Loop()
{
  if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      float kedep = 1000.0*edep;
      int ibin = int(floor((kedep-_emin)/_ebin));
      if(ibin>=0 && ibin<(int)_tdevec.size()){
	_evec[ibin]->Fill(kedep);
	_tdevec[ibin]->Fill(time_tcal-time_thv,mcshlen);
	_tdpevec[ibin]->Fill(time_tcal-time_thv,mcshlen);
      }
   }

   _tdecan = new TCanvas("tdecan","tdecan",700,700);
   unsigned nx = (unsigned)ceil(sqrt(_tdevec.size()));
   unsigned ny = ceil(_tdevec.size()/nx);
   TF1* line = new TF1("line","[0]+[1]*x",-10.0,10.0);
   line->SetParameter(0,0.0);
   line->FixParameter(0,0.0);
   line->SetParameter(1,100.0);
   _tdecan->Divide(nx,ny);
   for(unsigned iplot=0;iplot<_tdevec.size();++iplot){
     _tdecan->cd(iplot+1);
     _tdevec[iplot]->Draw("colorz");
     if(_tdpevec[iplot]->GetEntries() > 20){
       TFitResultPtr fitr = _tdpevec[iplot]->Fit(line,"QS","sames",-3.0,3.0);
       _smean.push_back(fitr->Parameter(1));
       _serr.push_back(fitr->ParError(1));
       _emean.push_back(_evec[iplot]->GetMean());
       _eerr.push_back(_evec[iplot]->GetRMS());
     }
   }

   _edcan = new TCanvas("edcan","edcan",700,700);
   _edcan->Divide(nx,ny);
   for(unsigned iplot=0;iplot<_evec.size();++iplot){
     _edcan->cd(iplot+1);
     _evec[iplot]->Draw();
   }
   
   _fits = new TCanvas("fits","fits",700,700);
   _fits->Divide(1,2);
   _fits->cd(1);
   TGraphErrors* slopes = new TGraphErrors(_emean.size(),_emean.data(),_smean.data(),_eerr.data(),_serr.data());
   slopes->SetName("TDEffV");
   slopes->SetTitle("(Half) TD Effective Propagation Velocity vs EDep;EDep (KeV);TD slope (mm/ns)");
   slopes->Draw("ALP");

   // now measure resolution
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      float kedep = 1000.0*edep;
      int ibin = min((int)_resvec.size()-1,max(0,int(floor((kedep-_emin)/_ebin))));
      float slope =_smean[ibin];
      if(ibin < (int)_resvec.size()-1)
	slope += (kedep-(_emin+ibin*_ebin))*(_smean[ibin+1]-_smean[ibin]);
      float tdlen = (time_tcal-time_thv)*slope;
      _resvec[ibin]->Fill(fabs(tdlen),tdlen-mcshlen);
      if(mcproc==56)_ceres->Fill(tdlen-mcshlen);
      if(mcpdg==2212)_pres->Fill(tdlen-mcshlen);
   }

   _rescan = new TCanvas("rescan","rescan",700,700);
   _rescan->Divide(nx,ny);
   for(unsigned iplot=0;iplot<_resvec.size();++iplot){
     _rescan->cd(iplot+1);
     _resvec[iplot]->Draw("colorz");
   }

   _resfitcan = new TCanvas("resfitcan","resfitcan",700,700);
   _resfitcan->Divide(nx,ny);
   TF1* rfit = new TF1("RFit",ALine,0.0,500.0,3);
   std::vector<float> p0;
   std::vector<float> p0err, p1err, p2err;
   for(unsigned iplot=0;iplot<_resvec.size();++iplot){
     _resfitcan->cd(iplot+1);
     _resvec[iplot]->FitSlicesY(0,0,-1,20);
     std::string fname(_resvec[iplot]->GetName());
     fname += "_2";
     TH1D *res_2 = (TH1D*)gDirectory->Get(fname.c_str());
     std::string ftitle(_resvec[iplot]->GetTitle());
     ftitle += "Fit Sigma";
     res_2->SetTitle(ftitle.c_str());
     res_2->SetLineColor(kRed);
     res_2->SetStats(1);
     res_2->SetMarkerStyle(1);
     rfit->SetParameters(_central,res_2->GetBinContent(1),0.05);
     rfit->FixParameter(0,_central);
     TFitResultPtr fitr = res_2->Fit(rfit,"SQ");
     p0.push_back(fitr->Parameter(0));
     p0err.push_back(fitr->ParError(0));
     _cres.push_back(fitr->Parameter(1));
     _creserr.push_back(fitr->ParError(1));
     _sres.push_back(fitr->Parameter(2));
     _sreserr.push_back(fitr->ParError(2));
   }
   TF1* rfun = new TF1("rfun","[0]+[1]*exp(-x/[2])",0.0,10.0);
   rfun->SetParameters(18,45.0,1.8);
   TF1* rsfun = new TF1("rsfun",DLine,0.0,10.0,4);
   rsfun->SetParameters(3.0,0.09,0.01,-0.02);
   _cresg = new TGraphErrors(_emean.size(),_emean.data(),_cres.data(),_eerr.data(),_creserr.data());
   _cresg->SetName("CentralRes");
   _cresg->SetTitle("TD Central Res vs EDep;EDep (KeV);Resolution (mm)");
   _sresg = new TGraphErrors(_emean.size(),_emean.data(),_sres.data(),_eerr.data(),_sreserr.data());
   _sresg->SetName("ResSlope");
   _sresg->SetTitle("TD Res Slope vs EDep;EDep (KeV);Resolution Slope");
   _rfits = new TCanvas("rfits","res fits", 700, 700);
   _rfits->Divide(2,2);
   _rfits->cd(1);
   _cresg->Fit(rfun,"Q","",0.5,6.5);
   _cresg->Draw("ALP");
   _rfits->cd(2);
   _sresg->Fit(rsfun,"Q","",0.5,6.5);
   _sresg->Draw("ALP");
   _rfits->cd(3);
   _ceres->SetLineColor(kRed);
   _ceres->Fit("gaus","Q");
   _rfits->cd(4);
   _pres->SetLineColor(kBlue);
   _pres->Fit("gaus","Q","same");
   TLegend* rleg = new TLegend(0.7,0.7,0.9,0.9);
   rleg->AddEntry(_ceres,"Ce","L");
   rleg->AddEntry(_pres,"Proton","L");

   // now pulls
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      float kedep = 1000.0*edep;
      int ibin = min((int)_resvec.size()-1,max(0,int(floor((kedep-_emin)/_ebin))));
      float slope =_smean[ibin];
      if(ibin < (int)_resvec.size()-1)
	slope += (kedep-(_emin+ibin*_ebin))*(_smean[ibin+1]-_smean[ibin]);
      float tdlen = (time_tcal-time_thv)*slope;
      float wlen = fabs(tdlen);
//      float tdres = _cres[ibin];
//      if(wlen > _central) tdres += _sres[ibin]*(wlen-_central);
      float tdres = rfun->Eval(kedep);
      if(wlen > _central) tdres += rsfun->Eval(kedep)*(wlen-_central);
      _pulledep->Fill(kedep,(tdlen-mcshlen)/tdres);
      _pullwlen->Fill(fabs(tdlen),(tdlen-mcshlen)/tdres);
      _pullsfrac->Fill(fabs(tdlen)/slen,(tdlen-mcshlen)/tdres);
      _pulltedep->Fill(1000.0*mcedep,(tdlen-mcshlen)/tdres);
      _pulltwlen->Fill(fabs(mcshlen),(tdlen-mcshlen)/tdres);
      _pulltsfrac->Fill(fabs(mcshlen)/slen,(tdlen-mcshlen)/tdres);
      if(mcproc==56)_cepull->Fill((tdlen-mcshlen)/tdres);
      if(mcpdg==2212)_ppull->Fill((tdlen-mcshlen)/tdres);
   }
   _pullcan = new TCanvas("pullcan","pullcan",1100,700);
   _pullcan->Divide(3,2);
   _pullcan->cd(1);
   _pulledep->Draw("colorz");
   _pullcan->cd(2);
   _pullwlen->Draw("colorz");
   _pullcan->cd(3);
   _pullsfrac->Draw("colorz");
   _pullcan->cd(4);
   _pulltedep->Draw("colorz");
   _pullcan->cd(5);
   _pulltwlen->Draw("colorz");
   _pullcan->cd(6);
   _pulltsfrac->Draw("colorz");

   float pullmin(0.2);
   float pullmax(2.0);
   _pulledep->FitSlicesY(0,0,-1,20);
   std::string pname(_pulledep->GetName());
   pname += "_2";
   TH1D *pulledep_2 = (TH1D*)gDirectory->Get(pname.c_str());
   pulledep_2->SetMinimum(pullmin);
   pulledep_2->SetMaximum(pullmax);
   pulledep_2->SetStats(0);
   pulledep_2->SetTitle("TD Pull Sigma vs EDep;EDep (KeV);Pull Sigma");
   pulledep_2->SetLineColor(kGreen);
   pulledep_2->SetMarkerColor(kGreen);
   pulledep_2->SetMarkerStyle(4);
   //
   _pulltedep->FitSlicesY(0,0,-1,20);
   std::string tpname(_pulltedep->GetName());
   tpname += "_2";
   TH1D *pulltedep_2 = (TH1D*)gDirectory->Get(tpname.c_str());
   pulltedep_2->SetMinimum(pullmin);
   pulltedep_2->SetMaximum(pullmax);
   pulltedep_2->SetStats(0);
   pulltedep_2->SetTitle("TD Pull Sigma vs EDep;EDep (KeV);Pull Sigma");
   pulltedep_2->SetLineColor(kBlack);
   pulltedep_2->SetMarkerColor(kBlack);
   pulltedep_2->SetMarkerStyle(5);
   //
   _pullwlen->FitSlicesY(0,0,-1,20);
   std::string wname(_pullwlen->GetName());
   wname += "_2";
   TH1D *pullwlen_2 = (TH1D*)gDirectory->Get(wname.c_str());
   pullwlen_2->SetMinimum(pullmin);
   pullwlen_2->SetMaximum(pullmax);
   pullwlen_2->SetStats(0);
   pullwlen_2->SetTitle("TD Pull Sigma vs Wire Length;Wire Length (mm);Pull Sigma");
   pullwlen_2->SetLineColor(kGreen);
   pullwlen_2->SetMarkerColor(kGreen);
   pullwlen_2->SetMarkerStyle(4);
   //
   _pulltwlen->FitSlicesY(0,0,-1,20);
   std::string twname(_pulltwlen->GetName());
   twname += "_2";
   TH1D *pulltwlen_2 = (TH1D*)gDirectory->Get(twname.c_str());
   pulltwlen_2->SetMinimum(pullmin);
   pulltwlen_2->SetMaximum(pullmax);
   pulltwlen_2->SetStats(0);
   pulltwlen_2->SetTitle("TD Pull Sigma vs Wire Length;Wire Length (mm);Pull Sigma");
   pulltwlen_2->SetLineColor(kBlack);
   pulltwlen_2->SetMarkerColor(kBlack);
   pulltwlen_2->SetMarkerStyle(5);
//
   _pullsfrac->FitSlicesY(0,0,-1,20);
   std::string sname(_pullsfrac->GetName());
   sname += "_2";
   TH1D *pullsfrac_2 = (TH1D*)gDirectory->Get(sname.c_str());
   pullsfrac_2->SetMinimum(pullmin);
   pullsfrac_2->SetMaximum(pullmax);
   pullsfrac_2->SetStats(0);
   pullsfrac_2->SetTitle("TD Pull Sigma vs Wire Fraction;Wire Fraction;Pull Sigma");
   pullsfrac_2->SetLineColor(kGreen);
   pullsfrac_2->SetMarkerColor(kGreen);
   pullsfrac_2->SetMarkerStyle(4);
   //
   _pulltsfrac->FitSlicesY(0,0,-1,20);
   std::string tsname(_pulltsfrac->GetName());
   tsname += "_2";
   TH1D *pulltsfrac_2 = (TH1D*)gDirectory->Get(tsname.c_str());
   pulltsfrac_2->SetMinimum(pullmin);
   pulltsfrac_2->SetMaximum(pullmax);
   pulltsfrac_2->SetStats(0);
   pulltsfrac_2->SetTitle("TD Pull Sigma vs Wire Fraction;Wire Fraction;Pull Sigma");
   pulltsfrac_2->SetLineColor(kBlack);
   pulltsfrac_2->SetMarkerColor(kBlack);
   pulltsfrac_2->SetMarkerStyle(5);
  //
   TLegend* pleg = new TLegend(0.7,0.7,0.9,0.9);
   pleg->AddEntry(pulledep_2,"Reco","PL");
   pleg->AddEntry(pulltedep_2,"MC True","PL");
   _pullfitcan = new TCanvas("pullfitcan","pullfitcan",1100,700);
   _pullfitcan->Divide(3,2);
   _pullfitcan->cd(1);
   pulledep_2->Draw();
   pulltedep_2->Draw("same");
   pleg->Draw();
   _pullfitcan->cd(2);
   pullwlen_2->Draw();
   pulltwlen_2->Draw("same");
   _pullfitcan->cd(3);
   pullsfrac_2->Draw();
   pulltsfrac_2->Draw("same");
   _pullfitcan->cd(4);
   _cepull->SetLineColor(kRed);
   _cepull->Fit("gaus","Q");
   _pullfitcan->cd(5);
   _ppull->SetLineColor(kBlue);
   _ppull->Fit("gaus","Q","same");
}

