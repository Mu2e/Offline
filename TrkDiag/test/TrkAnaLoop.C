//
// TrkAnaLoop.C
// -- example of how to loop through a TrkAna tree in a ROOT macro
// -- creates an example momentum, resolution and efficiency plot
// -- this is a work in progress

#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TChain.h"

#include "TMVA/Config.h"
#include "TMVA/Reader.h"

struct EventInfoTree {
  Int_t eventid, runid, subrunid;
  Int_t nprotons;
};

struct FitTree {
  Int_t status, pdg, nhits, ndof, nactive, ndouble, ndactive, nnullambig, nmat, nmatactive, nbend;
  Float_t t0, t0err, chisq, con, radlen, firstflt, lastflt, startvalid, endvalid,
    trkqual, mom, momerr, fltlen, d0, p0, om, z0, td, d0err, p0err, omerr, z0err, tderr;
};

struct MCEntTree {
  Float_t t0, mom, x, y, z, d0, p0, om, z0, td;
};

struct DeTrkqualTree{
  Double_t NActiveHits, ActiveHitFraction, Log10FitCon, MomError, T0Error, d0, MaxRadius, DoubleHitFraction, NullAmbigHitFraction, StrawHitFraction, trkqual;
};

void TrkAnaLoop(std::string filename, std::string treename) {

  // Get all the events
  double n_generated_ce_events = 0;
  double num_entries = 0;
  TChain* trkana = new TChain(treename.c_str());

  TFile* file = new TFile(filename.c_str(), "READ");
  if (!file->IsZombie()) {
    TH1F* hGenEventCount = (TH1F*) file->Get("genCountLogger/numEvents");
    if (!hGenEventCount) {
      std::cout << "Error: Couldn't get hGenEventCount" << std::endl;
      return;
    }
    n_generated_ce_events += hGenEventCount->GetBinContent(1);
    num_entries = hGenEventCount->Integral();
    trkana->Add(filename.c_str());
    file->Close();
  }
  std::cout << "CE Generated Events = " << n_generated_ce_events << std::endl;
  std::cout << "NEntries = " << num_entries << std::endl;


  // Create the histograms
  double min_mom = 101;
  double max_mom = 107;
  double bin_width = 0.05;
  int n_bins = (max_mom - min_mom) / bin_width;
  TH1F* hRecoCE = new TH1F("hRecoCE", "", n_bins,min_mom,max_mom);
  TH1F* hResCE = new TH1F("hResCE", "", 251,-4,4);

  TH1F* effPlot = new TH1F("CE Efficiency","",10,1,11);
  effPlot->GetXaxis()->SetBinLabel(1,"trkpt. && N_{digis} >= 10");
  effPlot->GetXaxis()->SetBinLabel(2,"pitch");
  effPlot->GetXaxis()->SetBinLabel(3,"d0");
  effPlot->GetXaxis()->SetBinLabel(4,"maxd");
  effPlot->GetXaxis()->SetBinLabel(5,"t0");
  effPlot->GetXaxis()->SetBinLabel(6,"trkqual");
  effPlot->GetXaxis()->SetBinLabel(7,"momentum");
  effPlot->SetTitle("CE Efficiency;Fraction of Events Passing Cut;Fraction of Generated Events");
  effPlot->SetMarkerSize(4);
  effPlot->GetXaxis()->SetTitleSize(.03);
  effPlot->GetYaxis()->SetTitleSize(.03);

  // Set the branch addresses for all the variables
  EventInfoTree evtinfo;
  FitTree de;
  MCEntTree demcent;
  Float_t pbi_evtwt;
  DeTrkqualTree detrkqual;
  trkana->SetBranchAddress("evtinfo.", &evtinfo);
  trkana->SetBranchAddress("de.", &de);
  trkana->SetBranchAddress("demcent", &demcent);
  trkana->SetBranchAddress("detrkqual", &detrkqual);
  trkana->SetBranchAddress("evtwt", &pbi_evtwt);

  // Book the TMVA reader
  TMVA::Reader* reader = new TMVA::Reader("V:Info");
  float nactive, frac_active, log_con, momerr, t0err, d0, max_radius, frac_dactive, frac_nullambig, frac_matactive;
  reader->AddVariable("detrkqual.NActiveHits", &nactive);
  reader->AddVariable("detrkqual.ActiveHitFraction", &frac_active);
  reader->AddVariable("detrkqual.Log10FitCon", &log_con);
  reader->AddVariable("detrkqual.MomError", &momerr);
  reader->AddVariable("detrkqual.T0Error", &t0err);
  reader->AddVariable("detrkqual.d0", &d0);
  reader->AddVariable("detrkqual.MaxRadius", &max_radius);
  reader->AddVariable("detrkqual.DoubleHitFraction", &frac_dactive);
  reader->AddVariable("detrkqual.NullAmbigHitFraction", &frac_nullambig);
  reader->AddVariable("detrkqual.StrawHitFraction", &frac_matactive);
  reader->BookMVA("MLP method", "dataset/TrkQualWeights/TMVAClassification_MLP.weights.xml");

  for (int i_entry = 0; i_entry < trkana->GetEntries(); ++i_entry) {

    trkana->GetEntry(i_entry);

    if (i_entry % 10000 == 0) {
      std::cout << i_entry << "/" << trkana->GetEntries() << std::endl;
    }

    double nominal_event_weight = pbi_evtwt;

    nactive = detrkqual.NActiveHits;
    frac_active = detrkqual.ActiveHitFraction;
    log_con = detrkqual.Log10FitCon;
    momerr = detrkqual.MomError;
    t0err = detrkqual.T0Error;
    d0 = detrkqual.d0;
    max_radius = detrkqual.MaxRadius;
    frac_dactive = detrkqual.DoubleHitFraction;
    frac_nullambig = detrkqual.NullAmbigHitFraction;
    frac_matactive = detrkqual.StrawHitFraction;

    double original_trk_qual = detrkqual.trkqual;
    double new_trk_qual = reader->EvaluateMVA("MLP method");

    
    effPlot->AddBinContent(1);
    if (de.td > 0.577350 && de.td < 1.000) {
      effPlot->AddBinContent(2);
      if (d0 > -80 && d0 < 105) {
        effPlot->AddBinContent(3);
        if (max_radius > 450 && max_radius < 680) {
          effPlot->AddBinContent(4);
          if (de.t0 > 700 && de.t0 < 1695) {
            effPlot->AddBinContent(5);
            if (new_trk_qual > 0.4) {
              effPlot->AddBinContent(6);
              if (de.mom > 103.85 && de.mom < 105.1) {
                effPlot->AddBinContent(7);
              }
            }
          }
        }
      }
    }
    //    std::cout << nactive << ", " << frac_active << ", " << log_con << ", " << momerr << ", " << t0err << ", " << d0 << ", " << max_radius << ", " << frac_dactive << ", " << frac_nullambig << ", " << frac_matactive << ", " << nominal_event_weight << ", " << new_trk_qual << ", " << original_trk_qual << std::endl;
    
    // The standard cuts
    if (de.t0>700 && de.t0<1695 && de.td > 0.577350 && de.td < 1.000 && d0>-80 && d0<105 && max_radius>450 && max_radius<680) {
      if (new_trk_qual > 0.4) {
	hRecoCE->Fill(de.mom, nominal_event_weight);
	hResCE->Fill(de.mom - demcent.mom, nominal_event_weight);
      }
    }
  }
 
  // Draw the histogram
  //hRecoCE->Draw("HIST E");
  effPlot->Scale(1/num_entries);
  effPlot->SetMaximum(1);
  effPlot->SetMinimum(0);
  //  effPlot->Draw("TEXT90"); 
}
