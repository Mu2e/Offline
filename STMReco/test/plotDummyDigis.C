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
  hRawWaveform->SetLineWidth(3);
  hRawWaveform->SetLineColor(kGreen);
  hRawWaveform->Draw("HIST");

  c1->cd(2);
  TH1F* hZSWaveform = (TH1F*) file->Get("plotZSWaveformDigis/evt0_waveform0");
  title = hZSWaveform->GetTitle();
  title += " ZS";
  hZSWaveform->SetTitle(title);
  hZSWaveform->SetLineWidth(3);
  hZSWaveform->SetLineColor(kBlue);
  hZSWaveform->Draw("HIST");

  c1->cd(3);
  TH1F* hPHWaveform = (TH1F*) file->Get("plotPHWaveformDigis/evt0_waveform0");
  title = hPHWaveform->GetTitle();
  title += " PH ";
  hPHWaveform->SetTitle(title);
  hPHWaveform->SetLineWidth(3);
  hPHWaveform->SetLineColor(kRed);
  hPHWaveform->Draw("HIST");

//New Canvas with all the histograms in the same axis
//No need to redefine the histograms from before
  TCanvas* c2 = new TCanvas();
  c2->SetTitle("STMWaveformDigis Super Imposed");
   
  hRawWaveform->SetLineColor(kGreen);
  hRawWaveform->SetLineWidth(3);

  hZSWaveform->SetLineColor(kBlue);
  hZSWaveform->SetLineWidth(3);

  hPHWaveform->SetLineColor(kRed);
  hPHWaveform->SetLineWidth(3);

  hRawWaveform->GetYaxis()->SetRangeUser(-1800,1800);
  hRawWaveform->Draw();
  hZSWaveform->Draw("SAME");
  hPHWaveform->Draw("SAME");

  auto legend = new TLegend();
  legend->AddEntry(hRawWaveform, "Raw","1");
  legend->AddEntry(hZSWaveform,"ZS","1");
  legend->AddEntry(hPHWaveform,"PH","1");
  //legend->Draw();

//New Canvas where the histogram plots the number of bins from the other two
  TCanvas* c3 = new TCanvas();
  c3->SetTitle("Number of bins from hPHWaveform and hZSWaveform");
  
  TH1F* hDigiSize = new TH1F("hDigiSize","Number of bins in histograms",1000,0,1000);
  hDigiSize->Fill(hPHWaveform->GetNbinsX());
  hDigiSize->Fill(hZSWaveform->GetNbinsX());
  hDigiSize->Draw();


 //Below is used to create a copy of the previous canvas
 // TCanvas* c4 = new TCanvas();
 // c4->DrawClonePad();

}

