//
//  Routine to calibrate the straw drift resolution
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
#include <ctime>

void DriftCalibRMS(TTree* ta,const char* dcut="",const char* ncut="") {
  if(ta==0){
    cout << "No tree specified" << endl;
    return;
  }

  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);

  double binsize = 0.025;
  double rdbinedges[2] = {0.0,2.5};
  unsigned finefactor(100);
  int nrdbins= (rdbinedges[1]-rdbinedges[0])/binsize;
  TH2D* sdr_vs_rd = new TH2D("sdr_vs_rd","Signed R_{drift}-R_{true} vs R_{drift};R_{drift} (mm);#Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1],120,-5.5,3.0);
  TProfile* sdr_vs_rdp = new TProfile("sdr_vs_rdp","Signed R_{drift}-R_{true} vs R_{drift};R_{drift} (mm);#Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1],-5.5,3.0,"s");
  TH2D* udr_vs_rd = new TH2D("udr_vs_rd","R_{drift}-abs(R_{true}) vs R_{drift};R_{drift} (mm);#Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1],120,-0.5,2.5);
  TProfile* udr_vs_rdp = new TProfile("udr_vs_rdp","R_{drift}-abs(R_{true}) vs R_{drift};R_{drift} (mm);#Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1],-0.5,2.5,"s");
  sdr_vs_rd->SetStats(0);
  sdr_vs_rdp->SetLineColor(kBlack);
  sdr_vs_rdp->SetMarkerColor(kBlack);
  udr_vs_rd->SetStats(0);
  udr_vs_rdp->SetLineColor(kRed);
  udr_vs_rdp->SetMarkerColor(kRed);

  TCut dsel(dcut);
  TCut nsel(ncut);
  TCut dhitsel =gfit+gmom+ghit+thit+dsel;
  TCut nhitsel =gfit+gmom+ghit+thit+nsel;
  ta->Project("sdr_vs_rd","copysign(detsh.rdrift,detsh.udoca*detshmc.doca)-detshmc.dist:detsh.rdrift",dhitsel);
  ta->Project("sdr_vs_rdp","copysign(detsh.rdrift,detsh.udoca*detshmc.doca)-detshmc.dist:detsh.rdrift",dhitsel);
  ta->Project("udr_vs_rd","detsh.rdrift-detshmc.dist:detsh.rdrift",nhitsel);
  ta->Project("udr_vs_rdp","detsh.rdrift-detshmc.dist:detsh.rdrift",nhitsel);

  TH1D* sdr_vs_rdp_mean = new TH1D("sdr_vs_rdp_mean","Mean #Delta R vs R_{drift};R_{drift} (mm);Mean #Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  TH1D* sdr_vs_rdp_sigma = new TH1D("sdr_vs_rdp_sigma","#sigma #Delta R vs R_{drift};R_{drift} (mm);#sigma #Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  sdr_vs_rdp_mean->SetLineColor(kBlack);
  sdr_vs_rdp_sigma->SetLineColor(kBlack);
  sdr_vs_rdp_sigma->SetStats(0);
  TH1D* udr_vs_rdp_mean = new TH1D("udr_vs_rdp_mean","Mean #Delta R vs R_{drift};R_{drift} (mm);Mean #Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  TH1D* udr_vs_rdp_sigma = new TH1D("udr_vs_rdp_sigma","#sigma #Delta R vs R_{drift};R_{drift} (mm);#sigma #Delta R (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  udr_vs_rdp_mean->SetStats(0);
  udr_vs_rdp_sigma->SetStats(0);
  udr_vs_rdp_mean->SetLineColor(kRed);
  udr_vs_rdp_sigma->SetLineColor(kRed);
  udr_vs_rdp_sigma->SetStats(0);
  auto lowbin = sdr_vs_rdp->FindBin(rdbinedges[0]+0.5*binsize);
  auto hibin = sdr_vs_rdp->FindBin(rdbinedges[1]-0.5*binsize);
  cout << "binedges [" << rdbinedges[0] <<"," << rdbinedges[1] << "] lowbin " << lowbin << " hibin " << hibin << endl;
  for(int ibin=1;ibin<=nrdbins;++ibin){
    sdr_vs_rdp_mean->SetBinContent(ibin,sdr_vs_rdp->GetBinContent(ibin));
    if(ibin < lowbin)
      sdr_vs_rdp_sigma->SetBinContent(ibin,sdr_vs_rdp->GetBinError(lowbin));
    else if(ibin > hibin)
      sdr_vs_rdp_sigma->SetBinContent(ibin,sdr_vs_rdp->GetBinError(hibin));
    else
      sdr_vs_rdp_sigma->SetBinContent(ibin,sdr_vs_rdp->GetBinError(ibin));
    udr_vs_rdp_mean->SetBinContent(ibin,udr_vs_rdp->GetBinContent(ibin));
    if(ibin < lowbin)
    udr_vs_rdp_sigma->SetBinContent(ibin,udr_vs_rdp->GetBinError(lowbin));
    else if(ibin > hibin)
    udr_vs_rdp_sigma->SetBinContent(ibin,udr_vs_rdp->GetBinError(hibin));
    else
    udr_vs_rdp_sigma->SetBinContent(ibin,udr_vs_rdp->GetBinError(ibin));
  }

  TLegend* cleg2 = new TLegend(0.3,0.1,0.7,0.4);
  cleg2->AddEntry(sdr_vs_rdp_sigma,"Signed #Delta R","PL");
  cleg2->AddEntry(udr_vs_rdp_sigma,"Unsigned #Delta R","PL");

  TCanvas* dcrmscan = new TCanvas("dcrmscan","dcrmscan",1000,1000);
  dcrmscan->Divide(2,2);
  dcrmscan->cd(1);
  gPad->SetLogz();
  sdr_vs_rd->Draw("colorz");
  sdr_vs_rdp->Draw("same");
  dcrmscan->cd(2);
  gPad->SetLogz();
  udr_vs_rd->Draw("colorz");
  udr_vs_rdp->Draw("same");
  dcrmscan->cd(3);
  sdr_vs_rdp_sigma->Draw();
  udr_vs_rdp_sigma->Draw("same");
  cleg2->Draw();
  //
  // dump the resolutions
  //
  string cfname;
  string cutstr = string(dcut)+string("_")+string(ncut);
  if(!cutstr.empty())
    cfname = string("DriftCalibRMS_") + cutstr + string(".txt");
  else
    cfname = string("DriftCalibRMS.txt");

  cout << "Saving calibration to " << cfname << endl;
  ofstream cfile(cfname.c_str(),ios::trunc);
  time_t now = time(0);
  char* dt = ctime(&now);
  cfile << "# The following was produced by DriftCalibRMS.C with drift, null hit selection " << cutstr << " on " << dt << endl;
  cfile << std::setw(4) << std::setprecision(3);
  cfile << "driftRMSBins : [ ";
  cfile << rdbinedges[0] << " , " << rdbinedges[1] << " ]" << endl;
  cfile << "unsignedDriftRMS : [ ";
  bool first = true;
  for(size_t ibin=0;ibin<nrdbins;++ibin){
    if(!first) cfile << " , ";
    first = false;
    cfile << udr_vs_rdp_sigma->GetBinContent(ibin+1);
  }
  cfile << " ]" << endl;
  cfile << "signedDriftRMS : [ ";
  first = true;
  for(size_t ibin=0;ibin<nrdbins;++ibin){
    if(!first) cfile << " , ";
    first = false;
    cfile << sdr_vs_rdp_sigma->GetBinContent(ibin+1);
 }
  cfile << " ]" << endl;
}

void DriftCalibRMSFile(const char* file,const char* dcut="",const char* ncut=""){
  TFile* tf = new TFile(file);
  TTree* ta = (TTree*)tf->Get("TAKK/trkana");
  DriftCalibRMS(ta,dcut,ncut);
}

void DriftCalibRMSChain(const char* files,const char* dcut="",const char* ncut=""){
  TChain* ta = new TChain("TAKK/trkana");
  FillChain(ta,files);
  DriftCalibRMS(ta,dcut,ncut);
}
