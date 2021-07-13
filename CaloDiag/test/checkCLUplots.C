//macro to perform basic checks on reco file
{
gStyle->SetOptStat(0);

TFile f("test.root");
f.cd("CaloExample");

TCanvas c1("c1","c1");
c1.Divide(3,2);

c1.cd(1);
Calo->Draw("cluEnergy>>h1");
h1->SetTitle("cluster energy");
h1->GetXaxis()->SetTitle("E (MeV)");
h1->SetLineWidth(2);
Calo->Draw("cluEnergy>>h2","cluConv==1","same");
h2->SetLineColor(2);
h2->SetLineWidth(2);
Calo->Draw("cluEnergy>>h2","cluConv==1","same");


c1.cd(4);
Calo->Draw("cluEnergyErr>>h3");
h3->SetTitle("cluster energy error");
h3->GetXaxis()->SetTitle("dE (MeV)");
h3->SetLineWidth(2);
Calo->Draw("cluEnergyErr>>h4","cluConv==1","same");
h4->SetLineColor(2);
h4->SetLineWidth(2);
Calo->Draw("cluEnergyErr>>h4","cluConv==1","same");


c1.cd(5);
Calo->Draw("cluTimeErr>>h5");
h5->SetTitle("cluster time error");
h5->GetXaxis()->SetTitle("dt (ns)");
h5->SetLineWidth(2);
Calo->Draw("cluTimeErr>>h6","cluConv==1","same");
h6->SetLineColor(2);
h6->SetLineWidth(2);
Calo->Draw("cluTimeErr>>h6","cluConv==1","same");

c1.cd(2);
Calo->Draw("cluTime>>h7");
h7->SetTitle("cluster time");
h7->GetXaxis()->SetTitle("t (ns)");
h7->SetLineWidth(2);
Calo->Draw("cluTime>>h8","cluConv==1","same");
h8->SetLineColor(2);
h8->SetLineWidth(2);
Calo->Draw("cluTime>>h8","cluConv==1","same");

c1.cd(3);
Calo->Draw("sqrt(cluCogX*cluCogX+cluCogY*cluCogY)>>h9");
h9->SetTitle("cluster radius");
h9->GetXaxis()->SetTitle("t (ns)");
h9->SetLineWidth(2);
Calo->Draw("sqrt(cluCogX*cluCogX+cluCogY*cluCogY)>>h10","cluConv==1","same");
h10->SetLineColor(2);
h10->SetLineWidth(2);
Calo->Draw("sqrt(cluCogX*cluCogX+cluCogY*cluCogY)>>h10","cluConv==1","same");


c1.cd(6);
Calo->Draw("cluNcrys>>h11");
h11->SetTitle("cluster #crystals");
h11->GetXaxis()->SetTitle("N crystals");
h11->SetLineWidth(2);
Calo->Draw("cluNcrys>>h12","cluConv==1","same");
h12->SetLineColor(2);
h12->SetLineWidth(2);
Calo->Draw("cluNcrys>>h12","cluConv==1","same");


c1.SaveAs("Cluster.gif");


}
