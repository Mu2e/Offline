// This script produces plots of Gaussian distributions of computed energy - Mcenergy


{
TCanvas* c1 = new TCanvas("c1", "ElectronFits", 400, 150);

c1->Divide(2,1);
c1->cd(1);

TFile f("electronFit25ns_8DoublePeak6Uniform.root");

treeDataElectronUniform = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
treeDataElectronUniform->Draw("(fittingFunction->GetParameter(1)/(6.8978069e6)-qMctrigenergy)*1000.0>>histElUniform(200,-1.0,1.0)");
histElUniform->Draw();
histElUniform->Fit("gaus");
gStyle->SetOptFit(1);
histElUniform->SetTitle("Difference of Reconstructed Energy and MCTrigenergy - 25 ns");
histElUniform->GetYaxis()->SetTitle("Frequency");
histElUniform->GetYaxis()->SetTitleOffset(1.4);
histElUniform->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MCTrigenergy [KeV]");

c1->cd(2);

TFile g("electronFit25ns_8DoublePeak6Gaussian.root");
treeDataElectronGaussian = (TTree*) gDirectory->Get("convolvedFitTree;1"); 
treeDataElectronGaussian->Draw("(fittingFunction->GetParameter(1)/(6.8978069e6)-qMctrigenergy)*1000.0>>histElGauss(200,-1.0,1.0)");
histElGauss->Draw();
histElGauss->Fit("gaus");
gStyle->SetOptFit(1);
histElGauss->SetTitle("Difference of Reconstructed Energy and MC Trigenergy - 25 ns");
histElGauss->GetYaxis()->SetTitle("Frequency");
histElGauss->GetYaxis()->SetTitleOffset(1.4);
histElGauss->GetXaxis()->SetTitle("Difference of Reconstructed Energy and MC Trigenergy [KeV]");



}