{
	TFile g("electronFit25ns_8UniformFirst.root");
	TTree *data = (TTree*) gDirectory->Get("convolvedFitTree");

	TF1 *func1;
	data->SetBranchAddress("fittingFunction",&func1);

	TFile f("electronFit25ns_8UniformSecond.root");
	TTree *data2 = (TTree*) gDirectory->Get("convolvedFitTree");
	
	TF1 *func2;
	data2->SetBranchAddress("fittingFunction",&func2);

	TH1F *diffChiHist = new TH1F("diffChiHist","diffChiHist",10000,-500,500);

	for (int i = 0; i < data->GetEntries(); ++i)
	{
		data->GetEntry(i);
		data2->GetEntry(i);
		if (func1->GetChisquare() != func2->GetChisquare())
		diffChiHist->Fill(func1->GetChisquare()-func2->GetChisquare());
	}

}