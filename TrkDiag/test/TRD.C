
#include "TrkDiag/test/FillChain.C+"
#include "TrkDiag/test/TrkRecoDiag.C+"

void TRD(const char* sigfiles,unsigned nsig, const char* bkgfiles) {
  TChain* mtce = new TChain("TrkRecoDiag/trdiag");
  FillChain(mtce,sigfiles);
  mtce->SetTitle("Sig ");
  mtce->SetLineColor(kRed);
  TChain* mtbk = new TChain("TrkRecoDiag/trdiag");
  FillChain(mtbk,bkgfiles);
  mtbk->SetTitle("Bkg ");
  mtbk->SetLineColor(kBlue);
  TrkRecoDiag trdce(mtce,nsig);
  trdce.Loop();
//  trdce.drawHistos();
  TrkRecoDiag trdbk(mtbk,mtbk->GetEntries());
  trdbk.Loop();
//  trdbk.drawHistos();
  TCanvas* acc = new TCanvas("acc","Acceptance",1000,1000);
  acc->Divide(1,2);
  acc->cd(1);
  trdce._acc->SetMaximum(1.0); // why is this necessary???
  trdce._acc->SetLabelSize(0.075);
  trdce._acc->SetMarkerSize(2.0);
  
  trdce._acc->Draw("PETEXT0");
  acc->cd(2);
  gPad->SetLogy();
  trdbk._acc->SetLabelSize(0.075);
  trdbk._acc->SetMarkerSize(2.0);
  trdbk._acc->Draw("PETEXT0");

// helix plots
  TLegend* hleg = new TLegend(0.6,0.7,0.8,0.9);
  hleg->AddEntry(trdce._hn,"Helix Fit","L");
  hleg->AddEntry(trdce._shn,"Chisq Fit","L");
  hleg->AddEntry(trdce._fhn,"Kalman Fit","L");
  hleg->AddEntry(trdce._phn,"Physics Sel","L");
//
  TCanvas* hel1 = new TCanvas("hel1","Helix 1",1000,1000);
  hel1->Divide(4,2);
  hel1->cd(1);
  trdce._hn->Draw();
  trdce._shn->Draw("same");
  trdce._fhn->Draw("same");
  trdce._phn->Draw("same");
  hel1->cd(2);
  trdce._hna->Draw();
  trdce._shna->Draw("same");
  trdce._fhna->Draw("same");
  trdce._phna->Draw("same");
  hleg->Draw();
  hel1->cd(3);
  trdce._hd0->Draw();
  trdce._shd0->Draw("same");
  trdce._fhd0->Draw("same");
  trdce._phd0->Draw("same");
  hel1->cd(4);
  trdce._hrmax->Draw();
  trdce._shrmax->Draw("same");
  trdce._fhrmax->Draw("same");
  trdce._phrmax->Draw("same");
//
  hel1->cd(5);
  gPad->SetLogy();
  trdbk._hn->Draw();
  trdbk._shn->Draw("same");
  trdbk._fhn->Draw("same");
  trdbk._phn->Draw("same");
  hel1->cd(6);
  gPad->SetLogy();
  trdbk._hna->Draw();
  trdbk._shna->Draw("same");
  trdbk._fhna->Draw("same");
  trdbk._phna->Draw("same");
  hel1->cd(7);
  gPad->SetLogy();
  trdbk._hd0->Draw();
  trdbk._shd0->Draw("same");
  trdbk._fhd0->Draw("same");
  trdbk._phd0->Draw("same");
  hel1->cd(8);
  gPad->SetLogy();
  trdbk._hrmax->Draw();
  trdbk._shrmax->Draw("same");
  trdbk._fhrmax->Draw("same");
  trdbk._phrmax->Draw("same");
//
  TCanvas* hel2 = new TCanvas("hel2","Helix 2",1000,1000);
  hel2->Divide(4,2);
  hel2->cd(1);
  trdce._hrad->Draw();
  trdce._shrad->Draw("same");
  trdce._fhrad->Draw("same");
  trdce._phrad->Draw("same");
  hel2->cd(2);
  trdce._hlam->Draw();
  trdce._shlam->Draw("same");
  trdce._fhlam->Draw("same");
  trdce._phlam->Draw("same");
  hel2->cd(3);
  trdce._hmom->Draw();
  trdce._shmom->Draw("same");
  trdce._fhmom->Draw("same");
  trdce._phmom->Draw("same");
//
  hel2->cd(5);
  gPad->SetLogy();
  trdbk._hrad->Draw();
  trdbk._shrad->Draw("same");
  trdbk._fhrad->Draw("same");
  trdbk._phrad->Draw("same");
  hel2->cd(6);
  gPad->SetLogy();
  trdbk._hlam->Draw();
  trdbk._shlam->Draw("same");
  trdbk._fhlam->Draw("same");
  trdbk._phlam->Draw("same");
  hel2->cd(7);
  gPad->SetLogy();
  trdbk._hmom->Draw();
  trdbk._shmom->Draw("same");
  trdbk._fhmom->Draw("same");
  trdbk._phmom->Draw("same");
// seed fit plots
  TLegend* sleg = new TLegend(0.6,0.7,0.8,0.9);
  sleg->AddEntry(trdce._ssna,"Chisq Fit","L");
  sleg->AddEntry(trdce._fsna,"Kalman Fit","L");
  sleg->AddEntry(trdce._psna,"Physics Sel","L");
  TCanvas* seed = new TCanvas("seed","Seed",1000,1000);
  seed->Divide(4,2);
  seed->cd(1);
  trdce._ssna->Draw();
  trdce._fsna->Draw("same");
  trdce._psna->Draw("same");
  sleg->Draw();
  seed->cd(2);
  trdce._ssmom->Draw();
  trdce._fsmom->Draw("same");
  trdce._psmom->Draw("same");
  seed->cd(3);
  trdce._ssmerr->Draw("same");
  trdce._fsmerr->Draw("same");
  trdce._psmerr->Draw("same");
  seed->cd(4);
  trdce._sschisq->Draw("same");
  trdce._fschisq->Draw("same");
  trdce._pschisq->Draw("same");
//
  seed->cd(5);
  gPad->SetLogy();
  trdbk._ssna->Draw();
  trdbk._fsna->Draw("same");
  trdbk._psna->Draw("same");
  seed->cd(6);
  gPad->SetLogy();
  trdbk._ssmom->Draw();
  trdbk._fsmom->Draw("same");
  trdbk._psmom->Draw("same");
  seed->cd(7);
  gPad->SetLogy();
  trdbk._ssmerr->Draw("same");
  trdbk._fsmerr->Draw("same");
  trdbk._psmerr->Draw("same");
  //
  seed->cd(8);
  gPad->SetLogy();
  trdbk._sschisq->Draw("same");
  trdbk._fschisq->Draw("same");
  trdbk._pschisq->Draw("same");
}
