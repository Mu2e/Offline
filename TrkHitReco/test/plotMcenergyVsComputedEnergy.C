{
		
	TFile f("TrkChargeReco/test/FitDataElectron.root");
	TTree* electronFitData = (TTree*) gDirectory->Get("FitTree");
	electronFitData->Draw("computedEnergy:mcenergy*1000","computedEnergy>0&&computedEnergy<1023&&mcenergy*1000<40");
	htemp->SetTitle("Reconstructed Energy vs. MC Energy");
	htemp->GetYaxis()->SetTitle("Reconstructed Energy (Counts)");
	htemp->GetXaxis()->SetTitle("MC Energy (KeV)");
	gPad->Update();
}	
