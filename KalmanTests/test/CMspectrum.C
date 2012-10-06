double CMSpectrum(double ee) {
  double mal(25133);
  double mmu(105.654);
  double emu(105.194);
  double emue(104.973);
  double me(0.511);
  double a5(8.6434e-17);
  double a6(1.16874e-17);
  double a7(-1.87828e-19);
  double a8(9.16327e-20);
  double delta = emu - ee - ee*ee/(2*mal);
  return a5*pow(delta,5) + a6*pow(delta,6) + a7*pow(delta,7) + a8*pow(delta,8);
}

Double_t CMspectrum(Double_t *x, Double_t *par) {
  Double_t ee = x[0];
  return par[0]*CMSpectrum(ee);
}

void plotCMSpectrum(double emin=85.0, double emax=104.973) {
  TF1* cm = new TF1("cm",CMspectrum,85.,105.,1);
  cm->SetParameter(0,1.0);
  double smax = CMSpectrum(emin);
  TH1F* cmspec = new TH1F("cmspec","CM Spectrum",100,emin,emax);
  cmspec->SetMaximum(smax);
  TCanvas* ccan = new TCanvas("ccan","ccan",900,600);
  cmspec->Draw();
  cm->Draw("same");

}
