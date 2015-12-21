void TestTrackerChargeDepositPlot() {

  gROOT->ProcessLine(".L TrackerChargeDepositPlots.C+");
  gROOT->ProcessLine(".L ConvEMinusRun.h+");
  gROOT->ProcessLine(".L BeamFlashRun.h+");

  ConvEMinusRun* conv = new ConvEMinusRun("~/Mu2e/Offline/data/batch/cd3-sims/IPA-Cone_conversion/*.root");
  BeamFlashRun* flash = new BeamFlashRun("~/Mu2e/Offline/data/batch/cd3-sims/IPA-All_flash/*.root");

  TrackerChargeDepositPlots("flash", flash, conv);
}
