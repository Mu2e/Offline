//
// This example ROOT macro plots the steps of the Moving Window Deconvolution algorithm
// This ROOT macro runs on the output of the STMMovingWindowDeconvolution module
// with verbosity level set to 5
//
void plotMWD(std::string filename = "mwd.root") {

  TFile* file = new TFile(filename.c_str(), "READ");

  TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.8);
  leg->SetTextSize(0.04);
  leg->SetLineColor(kWhite);

  std::string evt_wvf = "evt1_wvf161";

  TH1F* hRawWaveform = (TH1F*) file->Get(("mwdLaBr/h_waveform_"+evt_wvf).c_str());
  hRawWaveform->SetStats(false);
  hRawWaveform->SetLineColor(kBlack);
  hRawWaveform->Draw("HIST");
  leg->AddEntry(hRawWaveform, "raw waveform", "l");

  TH1F* hDeconvolved = (TH1F*) file->Get(("mwdLaBr/h_deconvolved_"+evt_wvf).c_str());
  hDeconvolved->SetLineColor(kOrange);
  hDeconvolved->Draw("HIST SAME");
  leg->AddEntry(hDeconvolved, "deconvolution (tau parameter)", "l");

  TH1F* hDifferentiated = (TH1F*) file->Get(("mwdLaBr/h_differentiated_"+evt_wvf).c_str());
  hDifferentiated->SetLineColor(kSpring);
  hDifferentiated->Draw("HIST SAME");
  leg->AddEntry(hDifferentiated, "differentiation (M parameter)", "l");


  TH1F* hAveraged = (TH1F*) file->Get(("mwdLaBr/h_averaged_"+evt_wvf).c_str());
  hAveraged->SetLineColor(kRed);
  hAveraged->SetLineWidth(2);
  hAveraged->Draw("HIST SAME");
  leg->AddEntry(hAveraged, "average (L parameter)", "l");

  TH1F* hBaselineMean = (TH1F*) file->Get(("mwdLaBr/h_baseline_mean_"+evt_wvf).c_str());
  hBaselineMean->SetLineColor(kBlue);
  hBaselineMean->SetLineWidth(2);
  hBaselineMean->Draw("HIST SAME");
  leg->AddEntry(hBaselineMean, "baseline mean", "l");

  TH1F* hPeaks = (TH1F*) file->Get(("mwdLaBr/h_peaks_"+evt_wvf).c_str());
  hPeaks->SetLineColor(kBlue);
  hPeaks->SetLineWidth(2);
  hPeaks->Draw("HIST SAME");
  leg->AddEntry(hPeaks, "found peaks", "l");

  TH1F* hPeakThresh = (TH1F*) file->Get(("mwdLaBr/h_peak_threshold_"+evt_wvf).c_str());
  hPeakThresh->SetLineColor(kRed);
  hPeakThresh->SetLineWidth(2);
  hPeakThresh->SetLineStyle(kDashed);
  hPeakThresh->Draw("HIST SAME");
  leg->AddEntry(hPeakThresh, "peak threshold", "l");

  double min_y = hRawWaveform->GetMinimum();
  if (hDeconvolved->GetMinimum() < min_y) {
    min_y = hDeconvolved->GetMinimum();
  }
  if (hDifferentiated->GetMinimum() < min_y) {
    min_y = hDifferentiated->GetMinimum();
  }
  if (hAveraged->GetMinimum() < min_y) {
    min_y = hAveraged->GetMinimum();
  }
  if (hPeaks->GetMinimum() < min_y) {
    min_y = hPeaks->GetMinimum();
  }

  double max_y = hRawWaveform->GetMaximum();
  if (hDeconvolved->GetMaximum() > max_y) {
    max_y = hDeconvolved->GetMaximum();
  }
  if (hDifferentiated->GetMaximum() > max_y) {
    max_y = hDifferentiated->GetMaximum();
  }
  if (hAveraged->GetMaximum() > max_y) {
    max_y = hAveraged->GetMaximum();
  }
  if (hPeaks->GetMaximum() > max_y) {
    max_y = hPeaks->GetMaximum();
  }
  hRawWaveform->GetYaxis()->SetRangeUser(min_y*1.1,max_y*1.1);

  leg->Draw();
}
