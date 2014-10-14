// This script produces plots of Gaussian distributions of computed energy - Mcenergy


{
	TCanvas* c1 = new TCanvas("c1", "ElectronFits", 400, 300);
	TCanvas* c2 = new TCanvas("c2", "ProtonFits", 400, 300);

c1->Divide(2,2);
c1->cd(1);
	TFile h("electronFit25ns_8DoublePeak.root"); 

	treeDataElectron25ns = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	treeDataElectron25ns->Draw("(func->GetParameter(1)/(6.8978069e6)-qMctrigenergy)*1000.0>>histEl25ns(200,-1.0,1.0)");
	histEl25ns->Draw();
	histEl25ns->Fit("gaus");
	gStyle->SetOptFit(1);
	histEl25ns->SetTitle("Difference of Reconstructed Energy and MCTrigenergy - 25 ns");
	histEl25ns->GetYaxis()->SetTitle("Frequency");
	histEl25ns->GetYaxis()->SetTitleOffset(1.4);
	histEl25ns->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MCTrigenergy [KeV]");

c1->cd(2);
treeDataElectron25ns->Draw("(func->GetParameter(1)/(6.8978069e6)-qMcenergy)*1000.0>>histEl25nsOrig(200,-1.0,1.0)");
histEl25nsOrig->Draw();
histEl25nsOrig->Fit("gaus");
gStyle->SetOptFit(1);
histEl25nsOrig->SetTitle("Difference of Reconstructed Energy and MC Energy - 25 ns");
histEl25nsOrig->GetYaxis()->SetTitle("Frequency");
histEl25nsOrig->GetYaxis()->SetTitleOffset(1.4);
histEl25nsOrig->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MC Energy [KeV]");

c1->cd(3);
	TFile k("electronFit25ns_8.root"); 

	treeDataElectron25ns2 = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	treeDataElectron25ns2->Draw("(func->GetParameter(1)/(6.8978069e6)-qMctrigenergy)*1000.0>>histEl25ns2(200,-1.0,1.0)");
	histEl25ns2->Draw();
	histEl25ns2->Fit("gaus");
	gStyle->SetOptFit(1);
	histEl25ns2->SetTitle("Difference of Reconstructed Energy and MCTrigenergy W/O shifted time - 25 ns");
	histEl25ns->GetYaxis()->SetTitle("Frequency");
	histEl25ns->GetYaxis()->SetTitleOffset(1.4);
	histEl25ns->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MCTrigenergy [KeV]");

	c1->cd(4);
treeDataElectron25ns2->Draw("(func->GetParameter(1)/(6.8978069e6)-qMcenergy)*1000.0>>histEl25nsOrig2(200,-1.0,1.0)");
histEl25nsOrig2->Draw();
histEl25nsOrig2->Fit("gaus");
gStyle->SetOptFit(1);
histEl25nsOrig2->SetTitle("Difference of Reconstructed Energy and MC Energy W/O shifted time - 25 ns");
histEl25nsOrig2->GetYaxis()->SetTitle("Frequency");
histEl25nsOrig2->GetYaxis()->SetTitleOffset(1.4);
histEl25nsOrig2->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MC Energy [KeV]");






//c1->SaveAs("electronFitsVariableSamplesTest.pdf");


c2->cd();
	TFile c("protonFit25ns_8.root"); 

	treeDataProton25ns = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
	treeDataProton25ns->Draw("(func->GetParameter(1)/(6.8692148e6)-qMcenergy)*1000.0>>histPr25ns(200,-1.0,1.0)");
	histPr25ns->Draw();
	histPr25ns->Fit("gaus");
	histPr25ns->SetTitle("Difference of Reconstructed Energy and MC Energy - 25 ns");
	histPr25ns->GetYaxis()->SetTitle("Frequency");
	histPr25ns->GetYaxis()->SetTitleOffset(1.4);
	histPr25ns->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MC Energy [KeV]");







}