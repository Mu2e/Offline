// Usage:
//
// root -l
// .x plotMuonTransport.C ("histMuonTransportSingleStage.root", "histMuonTransportStage2.root")
//
// The input trees are produced by MuonTransportSingleStage.fcl
// and MuonTransportStage1.fcl, MuonTransportStage2.fcl sequence.
//
// Andrei Gaponenko, 2013

//================================================================
TTree* getTree(const std::string& filename, const std::string& treename) {
  TFile *fff = new TFile(filename.c_str());
  TObject *obj = fff->Get(treename.c_str());
  if(!obj) {
    std::cerr<<"NULL for object name="<<treename<<" in file "<<filename<<std::endl;
    return 0;
  }
  TTree* tr = (TTree*)(obj);
  return tr;
}

//================================================================
TH1 *makeHisto(TTree *nt, const std::string& var, const std::string& selection, int nbins, double xmin, double xmax) {
  TH1 *hh = new TH1D(var.c_str(), var.c_str(), nbins, xmin, xmax);
  nt->Project(var.c_str(), var.c_str(), selection.c_str());
  hh->SetDirectory(0);
  return hh;
}

//================================================================
void plotMuonTransport(std::string singeStageFile, std::string secondStageFile) {
  TTree *ss = getTree(singeStageFile, "vdDumper/nt");
  TTree *s2 = getTree(secondStageFile, "vdDumper/nt");

  int vdfirst = 0;
  int vdlast = 15;

  TH1 *hs = makeHisto(ss, "volumeCopy", "(pdgId==13)", 1+vdlast-vdfirst, vdfirst-0.5, vdlast+0.5);
  hs->GetXaxis()->SetTitle("vd number");
  hs->GetYaxis()->SetTitle("muon count");
  hs->SetMinimum(0.);
  hs->SetMarkerStyle(20);

  TH1 *h2 = makeHisto(s2, "volumeCopy", "(pdgId==13)", 1+vdlast-vdfirst, vdfirst-0.5, vdlast+0.5);
  h2->SetLineColor(kRed);

  TCanvas *cc = new TCanvas();
  cc->SetGrid();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  hs->Draw("PE");
  h2->Draw("sames");
}

//================================================================
