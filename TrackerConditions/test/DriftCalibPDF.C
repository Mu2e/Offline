//
//  Routine to calibrate the straw drift correction using PDFs
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

void DriftCalibPDF(TTree* ta,const char* hcut="") {
  if(ta==0){
    cout << "No tree specified" << endl;
    return;
  }

  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.2);
  double binsize = 0.025;
  double doff = 10*binsize;
  double rdbinedges[2] = {0.0,2.5};
  unsigned finefactor(100);
  int nrdbins= (rdbinedges[1]-rdbinedges[0])/binsize;
  TH1D* mcdoca = new TH1D("mcdoca","Wire Distance;DOCA (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  double rcbinedges[2] = {rdbinedges[0]-doff,rdbinedges[1]+doff};
  int nrcbins= (rcbinedges[1]-rcbinedges[0])/binsize;
  TH1D* rcluster = new TH1D("rcluster","Cluster Drift Distance;R_{cluster} (mm)",nrcbins,rcbinedges[0],rcbinedges[1]);
  TH1D* rcluster_fine = new TH1D("rcluster_fine","Cluster Drift Distance;R_{cluster} (mm)",finefactor*nrcbins,rcbinedges[0],rcbinedges[1]);
  TH1D* rdrift = new TH1D("rdrift","Drift Distance;R_{drift} (mm)",nrdbins,rdbinedges[0],rdbinedges[1]);
  mcdoca->SetStats(0);
  rcluster->SetStats(0);
  rdrift->SetStats(0);

  TCut hsel(hcut);
  TCut hitsel =gfit+gmom+ghit+thit+hsel;
  ta->Project("mcdoca","detshmc.dist",hitsel);
  ta->Project("rcluster","detsh.cdrift",hitsel);
  ta->Project("rcluster_fine","detsh.cdrift",hitsel);
  ta->Project("rdrift","detsh.rdrift",hitsel);
  cout << "mcdoca entries " << mcdoca->GetEntries() << " rcluster entries " << rcluster->GetEntries() << endl;
  // normalize
  rcluster->Scale(1.0/rcluster->Integral());
  mcdoca->Scale(1.0/mcdoca->Integral());
  rdrift->Scale(1.0/rdrift->Integral());
  rdrift->SetMarkerColor(kGreen);
  rdrift->SetLineColor(kGreen);

  // loop over the bins in the measured pdf, and extract the shift necessary to make that replicate the true pdf

  std::vector<double> offsets(nrcbins,0.0);
  unsigned jbin(1);
  double nraw = 0.0;
  double ncalib = mcdoca->GetBinContent(jbin);
  double cedge= mcdoca->GetBinLowEdge(jbin);
  double calibsum = ncalib;
  double rawsum = nraw;
  for(int ibin=1;ibin<=nrcbins;++ibin){
    double bincount = rcluster->GetBinContent(ibin);
    double bincent = rcluster->GetBinCenter(ibin);
    rawsum += bincount;
    offsets[ibin-1] = bincent - (cedge + binsize*((nraw+0.5*bincount)/ncalib));
    if(nraw + bincount < ncalib){
      nraw += bincount;
    } else {
      // update for next bin
      nraw += bincount - ncalib;
      jbin++;
      if(jbin > nrdbins || mcdoca->GetBinContent(jbin) <= 0)break;
      // the remaining content goes into the next step
      ncalib = mcdoca->GetBinContent(jbin);
      cedge = mcdoca->GetBinLowEdge(jbin);
      calibsum += ncalib;
    }
  }
  // check
  if(jbin < nrcbins && (abs(1.0-calibsum)>1e-6 || abs(1.0-rawsum)>1e-6) ){
    cout << "early exit, jbin " << jbin << " remaining mc fraction " << 1.0-calibsum << " rcluster fraction " << 1.0-rawsum << endl;
  }

  TH1D* offsets_hist = new TH1D("offsets","R Offset vs R_{cluster};R_{cluster} (mm)",nrcbins,rcbinedges[0],rcbinedges[1]);
  offsets_hist->SetStats(0);
  TH1D* corrected = new TH1D("corrected","Corrected R_{cluster};Corrected R_{cluster} (mm)",nrcbins,rcbinedges[0],rcbinedges[1]);
  corrected->SetStats(0);
  corrected->SetLineColor(kRed);
  corrected->SetMarkerColor(kRed);
  offsets_hist->SetStats(0);
  for(size_t ibin=0;ibin<offsets.size();++ibin){
    offsets_hist->SetBinContent(ibin+1,offsets[ibin]);
  }
  for(size_t ibin=1;ibin <= finefactor*nrcbins;++ibin){
    unsigned jbin = floor(float(ibin)/float(finefactor));
    corrected->Fill(rcluster_fine->GetBinCenter(ibin+1)-offsets[jbin],rcluster_fine->GetBinContent(ibin+1));
  }
  corrected->Scale(1.0/corrected->Integral());

  TLegend* cleg = new TLegend(0.2,0.3,0.6,0.6);
  cleg->AddEntry(mcdoca,"MC DOCA","PL");
  cleg->AddEntry(rdrift,"Calibrated R_{drift}","PL");
  cleg->AddEntry(corrected,"Corrected R_{cluster}","PL");
   TCanvas* dccan = new TCanvas("dccan","dccan",1000,1000);
  dccan->Divide(2,2);
  dccan->cd(1);
  rcluster->Draw();
  dccan->cd(2);
  offsets_hist->Draw();
  dccan->cd(3);
  mcdoca->Draw();
  corrected->Draw("histsame");
  rdrift->Draw("same");
  cleg->Draw();
  //
  // dump the offsets
  //

  string cfname;
  string cutstr(hcut);
  if(!cutstr.empty())
    cfname = string("DriftCalibPDF_") + cutstr + string(".txt");
  else
    cfname = string("DriftCalibPDF.txt");

  cout << "Saving calibration to " << cfname << endl;
  ofstream cfile(cfname.c_str(),ios::trunc);
  time_t now = time(0);
  char* dt = ctime(&now);
  cfile << "# The following was produced by DriftCalibPDF.C with hit selection " << cutstr << " on " << dt << endl;
  cfile << std::setw(4) << std::setprecision(3);
  cfile << "driftOffBins : [ ";
  cfile << rcbinedges[0] << " , " << rcbinedges[1] << " ]" << endl;
  cfile << "driftOffset : [ ";
  bool first(true);
  for(size_t ibin=0;ibin<nrcbins;++ibin){
    if(!first) cfile << " , ";
    first = false;
    cfile << offsets[ibin];
  }
  cfile << " ]" << endl;
  cfile.close();
}

void DriftCalibPDFFile(const char* file,const char* hcut="") {
  TFile* tf = new TFile(file);
  TTree* ta = (TTree*)tf->Get("TAKK/trkana");
  DriftCalibPDF(ta,hcut);
}

void DriftCalibPDFChain(const char* files,const char* hcut="") {
  TChain* ta = new TChain("TAKK/trkana");
  FillChain(ta,files);
  DriftCalibPDF(ta,hcut);
}
