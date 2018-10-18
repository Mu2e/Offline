#include "TGraph.h"
#include "TH1F.h"

void computeRejectionGraph(TH1F &electronHist, TH1F &protonHist, const int numberOfBins, TGraph &rejectionGraph)
{
  Double_t truncX[numberOfBins], truncY[numberOfBins];

  int protonSum = 0;
  int electronSum = 0;
  for (int i = 1; i <= numberOfBins; ++i)
  {
    // For some reason bin number starts with 1 for TH1F
    protonSum += protonHist.GetBinContent(i);
    electronSum += electronHist.GetBinContent(i);
    // acceptance rate of electrons
    truncX[i - 1] = electronSum / (double) electronHist.GetEntries();
    // 1 - rejection rate
    truncY[i - 1] = 1 - (protonSum / (double) protonHist.GetEntries());
  }

  rejectionGraph = TGraph(numberOfBins,truncX,truncY);
}

void plotRejectionRate()
{
	const int numberOfBins = 400;
	TFile f("FitDataProtonMultiplePeaks.root");
	TTree* protonMultiplePeaksTree = (TTree*) gDirectory->Get("FitTree");
        TH1F protonMultiplePeaksHist = TH1F("protonMultiplePeaksHist","protonMultiplePeaksHist",numberOfBins,0.0,0.2); 
 
	TFile g("FitDataElectronMultiplePeaks.root");
	TTree* electronMultiplePeaksTree = (TTree*) gDirectory->Get("FitTree");
	TH1F electronMultiplePeaksHist = TH1F("electronMultiplePeaksHist","electronMultiplePeaksHist",numberOfBins,0.0,0.2);

	TGraph rejectionMultiplePeaksGraph;
	computeRejectionGraph(electronMultiplePeaksHist, protonMultiplePeaksHist, numberOfBins, rejectionMultiplePeaksGraph);
	rejectionMultiplePeaksGraph.Draw("apl");

/**
	TFile h("FitDataProtonSum.root");
	TFile k("FitDataElectronSum.root");
	
	TH1F electronHist25ns8 = TH1F("electronHist25ns8","electronHist25ns8",numberOfBins,0.0,0.2);
**/	 

}

