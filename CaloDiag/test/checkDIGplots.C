//Macro to perform basic checks on digi file
{
  gStyle->SetOptStat(0);

  TFile f("dig.brownd.CeEndpointMix.MDC2020d.001203_00000000.art");


  TCanvas c1("c1","c1");
  c1.Divide(2,2);

  c1.cd(1);
  Events->Draw("mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->pos_.x(): mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->pos_.y()>>h1","","colz");
  h1->GetXaxis()->SetTitle("x");
  h1->GetYaxis()->SetTitle("y");
  h1->SetTitle("CaloShowerStep position");

  c1.cd(2);
  Events->Draw("mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->pos_.z()>>h2");
  h2->GetXaxis()->SetTitle("z");
  h2->GetYaxis()->SetTitle("");
  h2->SetTitle("CaloShowerStep position");

  c1.cd(3)
  Events->Draw("mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->volumeG4ID()>>h3");
  h3->GetXaxis()->SetTitle("volume ID");
  h3->GetYaxis()->SetTitle("");
  h3->SetTitle("CaloShowerStep volumeID");

  c1.cd(4);
  Events->Draw("mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->time()>>h4","mu2e::CaloShowerSteps_compressDigiMCs__MixPrimary.product()->time()<2000");
  h4->GetXaxis()->SetTitle("time (ns)");
  h4->GetYaxis()->SetTitle("");
  h4->SetTitle("CaloShowerStep time");
  c1.SaveAs("DIG_check1.gif");



  TCanvas c2("c2","c2");

  Events->Draw("mu2e::CaloShowerSims_compressDigiMCs__MixPrimary.product()->time()>>hh1")
  hh1->GetXaxis()->SetTitle("time (ns)");
  hh1->GetYaxis()->SetTitle("");
  hh1->SetTitle("CaloShowerSim  Time");
  c2.SaveAs("DIG_check2.gif");



  TCanvas c3("c3","c3");
  c3.Divide(2,2);

  c3.cd(1);
  Events->Draw("mu2e::CaloDigis_CaloDigiMaker__MixPrimary.product()->t0()>>hh2")
  hh2->GetXaxis()->SetTitle("time (ns)");
  hh2->GetYaxis()->SetTitle("");
  hh2->SetTitle("CaloDigi  Time");

  c3.cd(2);
  Events->Draw("mu2e::CaloDigis_CaloDigiMaker__MixPrimary.product()->SiPMID()>>hh21")
  hh21->GetXaxis()->SetTitle("ID");
  hh21->GetYaxis()->SetTitle("");
  hh21->SetTitle("CaloDigi SiPM ID");
  c3.SaveAs("DIG_check3.gif");
}
