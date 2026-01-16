//
//
void plotDummyDigis(std::string filename = "stmWaveformDigis.root") {

  TFile* file = new TFile(filename.c_str(), "READ");

  TCanvas* c1 = new TCanvas();
  c1->Divide(2, 2);

  c1->cd(1);
  TH1F* hRawWaveform = (TH1F*) file->Get("plotRawWaveformDigis/evt0_waveform0");
  TString title(hRawWaveform->GetTitle());
  title += " RAW";
  hRawWaveform->SetTitle(title);
  hRawWaveform->Draw("HIST");

  c1->cd(2);
  TH1F* hZSWaveform = (TH1F*) file->Get("plotZSWaveformDigis/evt0_waveform0");
  title = hZSWaveform->GetTitle();
  title += " ZS";
  hZSWaveform->SetTitle(title);
  hZSWaveform->Draw("HIST");

  c1->cd(3);
  TH1F* hMWDWaveform = (TH1F*) file->Get("plotMWDWaveformDigis/evt0_waveform0");
  title = hMWDWaveform->GetTitle();
  title += " MWD";
  hMWDWaveform->SetTitle(title);
  hMWDWaveform->Draw("HIST");
}
