//
//  Routine to calibrate the straw drift correction and resolution (T2D) WRT the single-cluster drift function
//
#include "TCut.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TStyle.h"
#include <fstream>
#include "KKCuts.hh"
#include "FillChain.C"

void DriftCalibProf(TTree* ta,const char* hcut="") {
  if(ta==0){
    cout << "No tree specified" << endl;
    return;
  }

  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  double maxuderr(0.25);

  double binedges[2] = {-0.25,2.75};
  double binsize = 0.025;
  int NBINS= (binedges[1]-binedges[0])/binsize;

  TH2D* srt_vs_rc = new TH2D("srt_vs_rc","Signed R_{true} vs R_{cluster};R_{cluster} (mm);R_{true} (mm)",NBINS,binedges[0],binedges[1],240,-2.5,2.5);
  TH2D* urt_vs_rc = new TH2D("urt_vs_rc","Unsigned R_{true} vs R_{cluster};R_{cluster} (mm);R_{true} (mm)",NBINS,binedges[0],binedges[1],120,0.0,2.5);
  TH2D* ud_vs_rc = new TH2D("ud_vs_rc","Unsigned UDOCA vs R_{cluster};R_{cluster} (mm);abs(UDOCA) (mm)",NBINS,binedges[0],binedges[1],120,0.0,3.5);
  TH1D* uderr = new TH1D("uderr","UDOCA Error;UDOCA Error (mm)",100,0.0,1.0);
  TProfile* srt_vs_rcp = new TProfile("srt_vs_rcp","Signed R_{true} vs R_{cluster};R_{cluster} (mm);R_{true} (mm)",NBINS,binedges[0],binedges[1],-2.5,2.5,"s");
  TProfile* urt_vs_rcp = new TProfile("urt_vs_rcp","Unsigned R_{true} vs R_{cluster};R_{cluster} (mm);abs(R_{true}) (mm)",NBINS,binedges[0],binedges[1],-1,2.5,"s");
  TProfile* ud_vs_rcp = new TProfile("ud_vs_rcp","Unsigned UDOCA vs R_{cluster};R_{cluster} (mm);abs(UDOCA) (mm)",NBINS,binedges[0],binedges[1],0.0,3.5,"s");
  TProfile* uderr_vs_rcp = new TProfile("uderr_vs_rcp","UDOCA error vs R_{cluster};R_{cluster} (mm);UDOCA error (mm)",NBINS,binedges[0],binedges[1],0,maxuderr);
  urt_vs_rc->SetStats(0);
  ud_vs_rc->SetStats(0);
  srt_vs_rc->SetStats(0);
  srt_vs_rcp->SetLineColor(kRed);
  srt_vs_rcp->SetMarkerColor(kRed);
  urt_vs_rcp->SetLineColor(kRed);
  urt_vs_rcp->SetMarkerColor(kRed);
  ud_vs_rcp->SetMarkerColor(kRed);

  TCut hsel(hcut);
  char cutstring[60];
  snprintf(cutstring,60,"sqrt(detsh.udocavar)<%4.3f",maxuderr);
  TCut udocasel(cutstring);
  TCut signedsel =gfit+gmom+ghit+thit+dhit+hsel;
  TCut unsignedsel =gfit+gmom+ghit+thit+hsel;
  ta->Project("srt_vs_rc","copysign(detshmc.dist,detsh.udoca*detshmc.doca):detsh.cdrift",signedsel);
  ta->Project("srt_vs_rcp","copysign(detshmc.dist,detsh.udoca*detshmc.doca):detsh.cdrift",signedsel);
  ta->Project("urt_vs_rc","detshmc.dist:detsh.cdrift",unsignedsel);
  ta->Project("urt_vs_rcp","detshmc.dist:detsh.cdrift",unsignedsel);
  ta->Project("ud_vs_rc","abs(detsh.udoca):detsh.cdrift",unsignedsel+udocasel);
  ta->Project("ud_vs_rcp","abs(detsh.udoca):detsh.cdrift",unsignedsel+udocasel);
  ta->Project("uderr","sqrt(detsh.udocavar)",unsignedsel);
  ta->Project("uderr_vs_rcp","sqrt(detsh.udocavar):detsh.cdrift",unsignedsel);

  TH1D* srt_vs_rcp_sigma = new TH1D("srt_vs_rcp_sigma","#sigma R vs R_{cluster};R_{cluster} (mm);#sigma R (mm)",srt_vs_rcp->GetXaxis()->GetNbins(),srt_vs_rcp->GetXaxis()->GetXmin(),srt_vs_rcp->GetXaxis()->GetXmax());
  TH1D* urt_vs_rcp_sigma = new TH1D("urt_vs_rcp_sigma","#sigma R vs R_{cluster};R_{cluster} (mm);#sigma R (mm)",srt_vs_rcp->GetXaxis()->GetNbins(),srt_vs_rcp->GetXaxis()->GetXmin(),srt_vs_rcp->GetXaxis()->GetXmax());
  TH1D* ud_vs_rcp_sigma = new TH1D("ud_vs_rcp_sigma","#sigma R vs R_{cluster};R_{cluster} (mm);#sigma R (mm)",ud_vs_rcp->GetXaxis()->GetNbins(),ud_vs_rcp->GetXaxis()->GetXmin(),ud_vs_rcp->GetXaxis()->GetXmax());
  TH1D* ud_vs_rcp_sigma_corr = new TH1D("ud_vs_rcp_sigma_corr","#sigma R vs R_{cluster};R_{cluster} (mm);#sigma R (mm)",ud_vs_rcp->GetXaxis()->GetNbins(),ud_vs_rcp->GetXaxis()->GetXmin(),ud_vs_rcp->GetXaxis()->GetXmax());
  TH1D* srt_vs_rcp_mean = new TH1D("srt_vs_rcp_mean","Mean R vs R_{cluster};R_{cluster} (mm);Mean R (mm)",srt_vs_rcp->GetXaxis()->GetNbins(),srt_vs_rcp->GetXaxis()->GetXmin(),srt_vs_rcp->GetXaxis()->GetXmax());
  TH1D* urt_vs_rcp_mean = new TH1D("urt_vs_rcp_mean","Mean R vs R_{cluster};R_{cluster} (mm);Mean R (mm)",urt_vs_rcp->GetXaxis()->GetNbins(),urt_vs_rcp->GetXaxis()->GetXmin(),urt_vs_rcp->GetXaxis()->GetXmax());
  TH1D* ud_vs_rcp_mean = new TH1D("ud_vs_rcp_mean","Mean r vs R_{cluster};R_{cluster} (mm);Mean R (mm)",ud_vs_rcp->GetXaxis()->GetNbins(),ud_vs_rcp->GetXaxis()->GetXmin(),ud_vs_rcp->GetXaxis()->GetXmax());
  urt_vs_rcp_mean->SetLineColor(kRed);
  urt_vs_rcp_mean->SetMarkerColor(kRed);
  urt_vs_rcp_mean->SetStats(0);
  srt_vs_rcp_sigma->SetLineColor(kBlue);
  srt_vs_rcp_sigma->SetMarkerColor(kBlue);
  srt_vs_rcp_sigma->SetStats(0);
  urt_vs_rcp_sigma->SetLineColor(kRed);
  urt_vs_rcp_sigma->SetMarkerColor(kRed);
  urt_vs_rcp_sigma->SetStats(0);

  ud_vs_rcp_mean->SetLineColor(kGreen);
  ud_vs_rcp_mean->SetMarkerColor(kGreen);
  ud_vs_rcp_mean->SetStats(0);
  ud_vs_rcp_sigma->SetLineColor(kGreen);
  ud_vs_rcp_sigma->SetMarkerColor(kGreen);
  ud_vs_rcp_sigma->SetStats(0);
  ud_vs_rcp_sigma_corr->SetLineColor(kCyan);
  ud_vs_rcp_sigma_corr->SetMarkerColor(kCyan);
  ud_vs_rcp_sigma_corr->SetStats(0);

//  ud_vs_rc->FitSlicesY(0,1,NBINS);
//  TH1D* ud_vs_rc_fsy_mean = (TH1D*)gDirectory->Get("ud_vs_rc_1");
//  TH1D* ud_vs_rc_fsy_sigma = (TH1D*)gDirectory->Get("ud_vs_rc_2");
//  ud_vs_rc_fsy_mean->SetMarkerColor(kBlack);
//  ud_vs_rc_fsy_sigma->SetMarkerColor(kBlack);
//  ud_vs_rc_fsy_mean->SetLineColor(kBlack);
//  ud_vs_rc_fsy_sigma->SetLineColor(kBlack);

  for(int ibin=1;ibin<=srt_vs_rcp->GetXaxis()->GetNbins();++ibin){
    urt_vs_rcp_mean->SetBinContent(ibin,urt_vs_rcp->GetBinContent(ibin));
    urt_vs_rcp_sigma->SetBinContent(ibin,urt_vs_rcp->GetBinError(ibin));
    srt_vs_rcp_mean->SetBinContent(ibin,srt_vs_rcp->GetBinContent(ibin));
    srt_vs_rcp_sigma->SetBinContent(ibin,srt_vs_rcp->GetBinError(ibin));
    ud_vs_rcp_mean->SetBinContent(ibin,ud_vs_rcp->GetBinContent(ibin));
    ud_vs_rcp_sigma->SetBinContent(ibin,ud_vs_rcp->GetBinError(ibin));
    double rerr = ud_vs_rcp->GetBinError(ibin);
    double uderr = uderr_vs_rcp->GetBinContent(ibin);
    double rerr_corr = sqrt(max(0.0,rerr*rerr-uderr*uderr));
    cout << "udoca error raw value " << rerr << " average error " << uderr << " corrected value " << rerr_corr << endl;
    ud_vs_rcp_sigma_corr->SetBinContent(ibin,rerr_corr);
  }
  TLegend* dcleg = new TLegend(0.1,0.5,0.5,0.9);
  dcleg->AddEntry(srt_vs_rcp_sigma,"Signed R_{true}","L");
  dcleg->AddEntry(urt_vs_rcp_sigma,"abs(R_{true})","L");
  dcleg->AddEntry(ud_vs_rcp_sigma,"raw abs(UDOCA)","L");
  dcleg->AddEntry(ud_vs_rcp_sigma_corr,"corrected abs(UDOCA)","L");
//  dcleg->AddEntry(ud_vs_rc_fsy_sigma,"abs(UDOCA) FitSlicesY","Pl");

  TCanvas* dccan = new TCanvas("dccan","dccan",1500,1000);
  dccan->Divide(3,2);
  dccan->cd(1);
  gPad->SetLogz();
  urt_vs_rc->Draw("colorz");
  urt_vs_rcp->Draw("same");
  dccan->cd(2);
  gPad->SetLogz();
  srt_vs_rc->Draw("colorz");
  srt_vs_rcp->Draw("same");
  dccan->cd(3);
  TLine* udocaselg = new TLine(maxuderr,0.0,maxuderr,uderr->GetMaximum());
  udocaselg->SetLineColor(kRed);
  uderr->Draw();
  udocaselg->Draw();
  dccan->cd(4);
  gPad->SetLogz();
  ud_vs_rc->Draw("colorz");
  ud_vs_rcp->Draw("same");
  dccan->cd(5);
  urt_vs_rcp_mean->Draw();
  srt_vs_rcp_mean->Draw("same");
  ud_vs_rcp_mean->Draw("same");
//  ud_vs_rc_fsy_mean->Draw("same");
  dcleg->Draw();
  dccan->cd(6);
  ud_vs_rcp_sigma->Draw();
  srt_vs_rcp_sigma->Draw("same");
  urt_vs_rcp_sigma->Draw("same");
  ud_vs_rcp_sigma_corr->Draw("same");
//  ud_vs_rc_fsy_sigma->Draw("same");

  // dump the data within the calibration range
  //
  string cfname;
  string cutstr(hcut);
  if(!cutstr.empty())
    cfname = string("DriftCalibProf_") + cutstr + string(".txt");
  else
    cfname = string("DriftCalibrof.txt");

  cout << "Saving calibration to " << cfname << endl;
//  double crange[2] = {0.0,2.5};
  ofstream cfile(cfname.c_str(),ios::trunc);
  cfile << std::setw(4) << std::setprecision(3);
  cfile << "driftOffBins : [ ";
//  cfile << crange[0] << " , " << crange[1] << " ]" << endl;
  cfile << binedges[0] << " , " << binedges[1] << " ]" << endl;
  cfile << "driftOffset : [ ";
  bool first(true);
  for(size_t ibin=0;ibin<NBINS;++ibin){
    double rbin = urt_vs_rcp_mean->GetBinCenter(ibin);
//    if(rbin > crange[0] && rbin < crange[1]) {
      double elow = urt_vs_rcp_mean->GetBinLowEdge(ibin);
      if(!first) cfile << " , ";
      first = false;
      cfile << rbin-urt_vs_rcp_mean->GetBinContent(ibin+1);
//    }
  }
  cfile << " ]" << endl;
  cfile.close();
}

void DriftCalibProfFile(const char* file,const char* hcut="") {
  TFile* tf = new TFile(file);
  TTree* ta = (TTree*)tf->Get("TAKK/trkana");
  DriftCalibProf(ta,hcut);
}

void DriftCalibProfChain(const char* files,const char* hcut="") {
  TChain* ta = new TChain("TAKK/trkana");
  FillChain(ta,files);
  DriftCalibProf(ta,hcut);
}
