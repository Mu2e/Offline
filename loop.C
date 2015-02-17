{
	TFile f("electronFit25ns_8DoublePeak6Uniform.root");
	TTree *treeDataUniform = (TTree*) gDirectory->Get("convolvedFitTree");

	Double_t scalingUniform;
	Double_t scalingGauss;
	int option;
	TF1 *funcUniform = new TF1();
	TF1 *funcGauss = new TF1();

	treeDataUniform->SetBranchAddress("fittingFunction",&funcUniform);
	treeDataUniform->SetBranchAddress("option",&option);

	TFile g("electronFit25ns_8DoublePeak6Gaussian.root");
	TTree *treeDataGauss = (TTree*) gDirectory->Get("convolvedFitTree");
	treeDataGauss->SetBranchAddress("fittingFunction",&funcGauss);

	Double_t graphUniform[10000];
	Double_t graphGauss[10000];

	for (int i = 0; i < 10000; ++i )
	{
		treeDataGauss->GetEntry(i);
		treeDataUniform->GetEntry(i);
		if (option == 4 || option == 2)
		{
			graphGauss[i] = funcGauss->GetParameter(1);
			graphUniform[i] = funcUniform->GetParameter(1);
		}
		else
		{
			graphUniform[i] = 0.0;
			graphGauss[i] = 0.0;
		}
	}
	TGraph *gr = new TGraph(10000,graphGauss,graphUniform);
}