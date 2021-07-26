// Macro to perform basic checks on dts file
{
  gStyle->SetOptStat(0);

  TFile f("dts.brownd.CeEndpoint.MDC2020d.001203_00000000.art");

  TCanvas c1("c1","c1");
  c1.Divide(2,2);

  c1.cd(1);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->pos_.x(): mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->pos_.y()>>h1","","colz");
  h1->GetXaxis()->SetTitle("x");
  h1->GetYaxis()->SetTitle("y");
  h1->SetTitle("CaloShowerStep position");

  c1.cd(2);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->pos_.z()>>h2");
  h2->GetXaxis()->SetTitle("z");
  h2->GetYaxis()->SetTitle("");
  h2->SetTitle("CaloShowerStep position");

  c1.cd(3)->SetLogy(1);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->energyDepG4()>>h5");
  h5->GetXaxis()->SetTitle("energy Dep (MeV)");
  h5->GetYaxis()->SetTitle("");
  h5->SetTitle("CaloShowerStep energy Dep");

  c1.cd(4);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->time()>>h4","mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->time()<2500");
  h4->GetXaxis()->SetTitle("time (ns)");
  h4->GetYaxis()->SetTitle("");
  h4->SetTitle("CaloShowerStep time");
  c1.SaveAs("DST_check1.gif");


  TCanvas c2("c2","c2");
  c2.Divide(2,2);

  c2.cd(1)->SetLogy(1);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->nCompress()>>h3");
  h3->GetXaxis()->SetTitle("nCompress");
  h3->GetYaxis()->SetTitle("");
  h3->SetTitle("CaloShowerStep number compressed");

  c2.cd(2)->SetLogy(1);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->momentumIn()>>h6");
  h6->GetXaxis()->SetTitle("p_{In} (MeV)");
  h6->GetYaxis()->SetTitle("");
  h6->SetTitle("CaloShowerStep momentum In");

  c2.cd(3);
  Events->Draw("mu2e::CaloShowerSteps_compressDetStepMCs__Primary.product()->volumeG4ID()>>h7");
  h7->GetXaxis()->SetTitle("vol ID");
  h7->GetYaxis()->SetTitle("");
  h7->SetTitle("CaloShowerStep volume ID");
  c2.SaveAs("DST_check2.gif");

}
