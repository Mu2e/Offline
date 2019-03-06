#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include <iostream>
#include <fstream>

// characterize the hit error as a function of the DOCA
void FitHitErrors(TTree* ta) {
// pure radial error
  TProfile* dc = new TProfile("dc","Unsigned Drift Calibration;DOCA (mm);Unsigned Drift Residual (mm)",50,0,2.5,-1.5,1.5,"s");
  TH2F* dch = new TH2F("dch","Unsigned Drift Calibration;DOCA (mm);Unsigned Drift Residual (mm)",50,0,2.5,100,-3.0,3.0);
// error including ambiguity assignment 
  TProfile* dcs = new TProfile("dcs","Signed Drift Calibration;DOCA (mm);Signed Drift Residual (mm)",50,0,2.5,-1.5,1.5,"s");
  TH2F* dchs = new TH2F("dchs","Signed Drift Calibration;DOCA (mm);Signed Drift Residual (mm)",50,0,2.5,100,-3.0,3.0);
  TCut goodtrack("de.status>0&&de.trkqual>0.2");
  TCut goodmctrack("demc.gen==2&&demcent.mom>100&&demcent.mom-demcxit.mom<1.5");
  TCut goodhit("detsh._active&&detshmc._rel==0");
  ta->Project("dc","detsh._rdrift-detshmc._dist:min(abs(detsh._doca),2.499)",goodtrack&&goodmctrack&&goodhit);
  ta->Project("dch","detsh._rdrift-detshmc._dist:min(abs(detsh._doca),2.499)",goodtrack&&goodmctrack&&goodhit);
  ta->Project("dcs","detsh._rdrift*detsh._ambig-detshmc._dist*detshmc._ambig:min(abs(detsh._doca),2.499)",goodtrack&&goodmctrack&&goodhit);
  ta->Project("dchs","detsh._rdrift*detsh._ambig-detshmc._dist*detshmc._ambig:min(abs(detsh._doca),2.499)",goodtrack&&goodmctrack&&goodhit);

  TCanvas* dcan =  new TCanvas("dcan","dcan",800,800);
  dcan->Divide(2,2);
  dcan->cd(1);
  dch->Draw("colorz");
  dc->Draw("same");
  dcan->cd(2);
  dc->Draw();
  dcan->cd(3);
  dchs->Draw("colorz");
  dcs->Draw("same");
  dcan->cd(4);
  dcs->Draw();

  unsigned nx = dc->GetNbinsX();
  std::vector<double> err(nx,0.0);
  std::vector<double> mean(nx,0.0);
  std::vector<double> errs(nx,0.0);
  std::vector<double> doca(nx,0.0);
 

  for(unsigned ibin=1;ibin<=nx;++ibin){
    err[ibin-1] = dc->GetBinError(ibin);
    mean[ibin-1] = dc->GetBinContent(ibin);
    errs[ibin-1] = dcs->GetBinError(ibin);
    doca[ibin-1] = dc->GetBinCenter(ibin);

  }
  std::filebuf fb;
  fb.open ("unsignederr.txt",std::ios::out);
  std::ostream os(&fb);
  for(auto errval : err){
    os << errval << ", ";
  }
  os << std::endl;
  fb.close();
  fb.open ("signederr.txt",std::ios::out);
  std::ostream os2(&fb);
  for(auto errval : errs){
    os2 << errval << ", ";
  }
  os2 << std::endl;
  fb.close();



  TGraph* ed = new TGraph(nx,doca.data(),err.data());
  ed->SetTitle("Average Hit Unsigned Radius RMS vs DOCA;DOCA (mm);R RMS (mm)");
  ed->SetMarkerStyle(20);
  TGraph* md = new TGraph(nx,doca.data(),mean.data());
  md->SetTitle("Average Hit Unsigned Radius vs DOCA;DOCA (mm);R RMS (mm)");
  md->SetMarkerStyle(20);
  TGraph* eds = new TGraph(nx,doca.data(),errs.data());
  eds->SetTitle("Average Hit Signed Radius RMS vs DOCA;DOCA (mm);R RMS (mm)");
  eds->SetMarkerStyle(20);

  TCanvas* fcan = new TCanvas("fcan","fcan",800,800);
  fcan->Divide(2,2);
  fcan->cd(1);
  md->Draw("AC*");
  fcan->cd(2);
  ed->Draw("AC*");
  fcan->cd(3);
  eds->Draw("AC*");


}
