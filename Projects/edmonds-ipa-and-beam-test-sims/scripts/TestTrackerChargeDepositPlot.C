#include "IPASims_LinkDef.h"

void TestTrackerChargeDepositPlot() {

  gROOT->ProcessLine(".L TrackerChargeDepositPlots.C+");
  gROOT->ProcessLine(".L BaseRun.h+");
  gROOT->ProcessLine(".L ConvEMinusRun.h+");
  gROOT->ProcessLine(".L BeamFlashRun.h+");
  gROOT->ProcessLine(".L ProtonRun.h+");

  ConvEMinusRun* conv = new ConvEMinusRun("~/Mu2e/Offline/data/batch/cd3-sims/IPA-Cone_conversion/*.root");
  BeamFlashRun* flash = new BeamFlashRun("~/Mu2e/Offline/data/batch/cd3-sims/IPA-All_flash/*.root");
  ProtonRun* proton   = new ProtonRun("~/Mu2e/Offline/data/batch/cd3-sims/IPA-Cone_proton/*.root");

  std::vector<BaseRun*> runs;

  runs.push_back(flash);
  TrackerChargeDepositPlots("flash", runs, conv);

  runs.push_back(proton);
  TrackerChargeDepositPlots("flash_proton", runs, conv);

  double integral = hDeadStrawCms_flash_proton->Integral() * (hDeadStrawCms_flash_proton->GetBinWidth(1) / (hDeadStrawCms_flash_proton->GetXaxis()->GetXmax() - hDeadStrawCms_flash_proton->GetXaxis()->GetXmin()));
  std::cout << "AE: " << integral << std::endl;

  hDeadStrawCms_flash_proton->Draw("HIST");
  TLine* line = new TLine(hDeadStrawCms_flash_proton->GetXaxis()->GetXmin(), integral, hDeadStrawCms_flash_proton->GetXaxis()->GetXmax(), integral);
  line->Draw("LSAME");
}
