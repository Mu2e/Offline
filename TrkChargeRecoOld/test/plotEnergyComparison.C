{
	TFile f("Mu2e3/test/FitDataElectronSum.root");
	TTree* electronFitData = (TTree*) gDirectory->Get("FitTree");
	electronFitData->Draw("computedEnergy>>histEl(400,0,10000)","computedEnergy>0&&computedEnergy<10000");
	histEl->SetLineColor(kBlue);

	TFile g("Mu2e3/test/FitDataProtonSum.root");
	TTree* protonFitData = (TTree*) gDirectory->Get("FitTree");
	protonFitData->Draw("computedEnergy>>histPr(400,0,10000)","computedEnergy>0&&computedEnergy<10000","same");
	histPr->SetLineColor(kRed);

	histEl->SetTitle("Reconstructed Energies");
	histEl->GetXaxis()->SetTitle("Reconstructed Energy (counts)");
	histEl->Draw();
	histPr->Draw("same");


	TLegend *legend = new TLegend(.6,.65,.88,.85);
	legend->AddEntry(histEl,"Electron","lpe");
	legend->AddEntry(histPr,"Proton","lpe");
	legend->Draw();
}
