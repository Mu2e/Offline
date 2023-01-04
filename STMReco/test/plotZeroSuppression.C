//
// This example ROOT macro plots an unsuppressed STMWaveformDigi (blue) and a handful of zero-suppressed digis on top (red).
// This ROOT macro runs on the output of Offline/STMReco/fcl/plotSTMWaveformDigis.fcl
//
void plotZeroSuppression(std::string filename = "stmWaveformDigis.root", int n_zp_digis_to_draw = 5) {

  TFile* file = new TFile(filename.c_str(), "READ");

  TH1F* hRawWaveform = (TH1F*) file->Get("plotHPGeWaveformDigis/evt1_waveform0");
  hRawWaveform->Draw("HIST");

  for (int i_zp_digi = 0; i_zp_digi < n_zp_digis_to_draw; ++i_zp_digi) {
    TH1F* hZPWaveform = (TH1F*) file->Get(Form("plotHPGeWaveformDigisZP/evt1_waveform%d", i_zp_digi));
    hZPWaveform->SetLineColor(kRed);
    hZPWaveform->Draw("HIST SAME");
  }
}
