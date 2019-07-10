void TrkAnaTutEx02() {
  std::string filename = "trkana-ex02.root";
  std::string treename = "TrkAnaEx02/trkana";

  TFile* file = new TFile(filename.c_str(), "READ");
  TTree* trkana = (TTree*) file->Get(treename.c_str());

  double min_mom = 101;
  double max_mom = 106;
  double mom_width = 0.5;
  int n_mom_bins = (max_mom - min_mom) / mom_width;
  std::string histname = "hRecoMom";
  TH1F* hRecoMom = new TH1F(histname.c_str(), "Reconstructed Momentum", n_mom_bins,min_mom,max_mom);

  std::string drawcmd = "trkent.mom>>" + histname;
  trkana->Draw(drawcmd.c_str(), "");
}
