// Optimimze Time Component of Fit
{
	TFile h("electronFit25ns_8DoublePeak6.root"); 
	treeDataElectron25ns = (TTree*) gDirectory->Get("convolvedFitTree;1"); 

	TF1 *func;
	Double_t measurementTimes[8];
	Double_t adc[8];
	TF1 *func7;
	TF1 *func8;

	treeDataElectron25ns->SetBranchAddress("qAdc",adc);
	treeDataElectron25ns->SetBranchAddress("func",&func);
	treeDataElectron25ns->SetBranchAddress("func7",&func7);
	treeDataElectron25ns->SetBranchAddress("func8",&func8);
	treeDataElectron25ns->SetBranchAddress("qMeasurementTimes",measurementTimes);

	Double_t qAdc[8];
	Double_t qMeasurementTimes[8];
	TF1 *qFunc;
	Int_t locMax;
	float probChi7, probChi8;

	TFile *newfile = new TFile("temp.root","recreate");
	TTree *temp = new TTree("tempTree","tempTree");
	temp->Branch("qAdc",qAdc,"qAdc[8]/D");
	temp->Branch("qMeasurementTimes",qMeasurementTimes,"qMeasurementTimes[8]/D");
	temp->Branch("func",&qFunc);
	temp->Branch("locMax",&locMax,"locaMax/I");
	temp->Branch("probChi7",&probChi7);
	temp->Branch("probChi8",&probChi8);


	Double_t locMaxArray[1e4];
	Double_t shiftTime[1e4];


	for (int i = 0; i < 1.0e4; ++i)
	{
		treeDataElectron25ns->GetEntry(i);
		temp->GetEntry(i);
		for (int j= 0; j < 8; ++j)
		{
		qAdc[j] = adc[j];
		qMeasurementTimes[j] = measurementTimes[j];
		}
		locMax = TMath::LocMax(8,qAdc);
		locMaxArray[i] = locMax;
		qFunc = func;
		shiftTime[i] = func->GetParameter(0);
		probChi8 = TMath::Prob(func8->GetChisquare(),3);
		probChi7 = TMath::Prob(func7->GetChisquare(),4);
		temp->Fill();
	}
	temp->AutoSave();
}