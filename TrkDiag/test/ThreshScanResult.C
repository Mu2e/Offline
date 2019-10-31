#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
void DrawThresh() {
  std::vector<double> thresh {5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  std::vector<double> totacc {0.2053,0.2203, 0.2233, 0.2176, 0.2098, 0.1987, 0.1607};
  std::vector<double> fitacc {0.7581,0.8242, 0.8341, 0.8280, 0.8401, 0.8253, 0.7736};
  std::vector<double> qualacc {0.9078,0.9064, 0.9038, 0.8933, 0.8667, 0.8344, 0.7363};
  std::vector<double> sigma { 0.141627, 0.148026, 0.14388, 0.153407, 0.153914, 0.171091, 0.199092};
  std::vector<double> sigerr { 0.0046428, 0.00374016, 0.00455301, 0.00523272, 0.0055792, 0.00494009, 0.0055657};
  std::vector<double> ppos{ 7.94547, 6.60994, 10.0901, 11.8033, 18.0821, 14.3907, 19.0542};
  std::vector<double> pposerr{ 2.11332, 1.62303, 3.30519, 5.43639, 10.2412, 8.49268, 16.1213};

  std::vector<double> ofrac{ 0.00100597, 0.000403619, 0.000408514, 0.000426369, 0.000554674, 0.000102494, 0.000499756};

  TGraph* tacc = new TGraph(thresh.size(),thresh.data(),totacc.data());
  TGraph* facc = new TGraph(thresh.size(),thresh.data(),fitacc.data());
  TGraph* qacc = new TGraph(thresh.size(),thresh.data(),qualacc.data());

  tacc->SetLineColor(kRed);
  tacc->SetMarkerColor(kRed);
  tacc->SetMarkerStyle(20);
  facc->SetLineColor(kBlue);
  facc->SetMarkerColor(kBlue);
  facc->SetMarkerStyle(21);
  qacc->SetLineColor(kGreen);
  qacc->SetMarkerColor(kGreen);
  qacc->SetMarkerStyle(22);

  TGraphErrors* sigg = new TGraphErrors(thresh.size(),thresh.data(),sigma.data(),0,sigerr.data());
  sigg->SetLineColor(kBlack);
  sigg->SetMarkerStyle(20);
  sigg->SetMarkerColor(kBlack);

  TGraphErrors* pposg = new TGraphErrors(thresh.size(),thresh.data(),ppos.data(),0,pposerr.data());
  pposg->SetLineColor(kBlack);
  pposg->SetMarkerStyle(21);
  pposg->SetMarkerColor(kBlack);

  TGraph* ofracg = new TGraph(thresh.size(),thresh.data(),ofrac.data());
  ofracg->SetLineColor(kBlack);
  ofracg->SetMarkerColor(kBlack);
  ofracg->SetMarkerStyle(22);


  TH1F* aacc = new TH1F("aacc","Absolute Track Acceptance vs Threshold;Thresold (mV);Acceptance",100,0.0,25.0);
  TH1F* racc = new TH1F("racc","Relative Track Acceptance vs Threshold;Thresold (mV);Acceptance",100,0.0,25.0);
  aacc->SetMinimum(0.15);
  aacc->SetMaximum(0.25);
  aacc->SetStats(0);
  racc->SetMinimum(0.7);
  racc->SetMaximum(1.0);
  racc->SetStats(0);
  TCanvas* acan = new TCanvas("acan","acan",800,800);
  acan->Divide(1,2);
  acan->cd(2);
  aacc->Draw();
  tacc->Draw("pLs");
  acan->cd(1);
  racc->Draw();
  facc->Draw("pLs");
  qacc->Draw("pLs");
  TLegend* aleg = new TLegend(0.6,0.7,0.9,0.9);
  aleg->AddEntry(tacc,"Total Acceptance","LP");
  aleg->AddEntry(facc,"Fit Acceptance","LP");
  aleg->AddEntry(qacc,"Quality Acceptance","LP");
  aleg->Draw();

  TH1F* sig = new TH1F("sig","Core Resolution Sigma vs Threshold;Thresold (mV);Sigma (MeV)",100,0.0,25.0);
  sig->SetMinimum(0.1);
  sig->SetMaximum(0.25);
  sig->SetStats(0);
  TH1F* pposh = new TH1F("pposh","Positive Tail Power vs Threshold;Thresold (mV);Power",100,0.0,25.0);
  pposh->SetMinimum(0.0);
  pposh->SetMaximum(25.0);
  pposh->SetStats(0);
  TH1F* ofrach = new TH1F("ofrach","Outlier Fraction vs Threshold;Thresold (mV);Fraction",100,0.0,25.0);
  ofrach->SetMinimum(0.0);
  ofrach->SetMaximum(0.0015);
  ofrach->SetStats(0);
  TCanvas* rcan = new TCanvas("rcan","rcan",800,800);
  rcan->Divide(2,2);
  rcan->cd(1);
  sig->Draw();
  sigg->Draw("pLs");
  rcan->cd(2);
  pposh->Draw();
  pposg->Draw("pLs");
  rcan->cd(3);
  ofrach->Draw();
  ofracg->Draw("pLs");


}
