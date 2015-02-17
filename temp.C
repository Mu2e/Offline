{

// Look at data with many more samples in PC_Waveform.root
//void fitData(TString dataName)
//{
	TString dataName = "SWF396_1;1";

	TFile f("PC_waveform.root");
	f.cd("makeSD");

	const double shapingTime = 25.0;
	cout << dataName << endl;
	TH1F *histData = (TH1F*) gDirectory->Get(dataName); 

	int option = 0;
	TF1 *fittingFunction = new TF1();
	
	std::vector<unsigned int> *adc = new std::vector<unsigned int>;
	float mcenergy;
	float mctrigenergy;

	const int numEntries = histData->GetEntries();
	Double_t qMeasurementTimes[numEntries];
	Double_t qAdc[numEntries];
	Double_t errorY[numEntries];
	Double_t errorX[numEntries];

	float qMcenergy;
	float qMctrigenergy;

	// Convert histogram to TGraph
	for (int i = 0; i < histData->GetEntries(); ++i)
	{
		qAdc[i] = histData->GetBinContent(i);
		qMeasurementTimes[i] = histData->GetBinCenter(i);
		errorY[i] = 0.75;
		errorX[i] = 0.0;
	}

	TGraphErrors *graph = new TGraphErrors(histData->GetEntries(),qMeasurementTimes,qAdc,errorX,errorY);


		std::vector<Float_t> tPeak;
		std::vector<Float_t> adcPeak;
		findPeaks(graph,tPeak,adcPeak,2.0);
			if (tPeak.size() == 0){
			cout << graph->GetX()[0] << endl;
			cout << "fail " << endl;}

		fittingFunction = new TF1();


		const double firstBin = graph->GetX()[0];
		const double lastBin = graph->GetX()[graph->GetN()-1];
		// If we have a dynamic pedestal
		if (tPeak[0] == 0.0)
		{
			if (tPeak.size() == 1)
			{
				option = 1;

				const double Q = TMath::Max(qAdc[0],0.0);
				Double_t newData[8];

				// Subtract dynamic pedestal from data
				for (int i = 0; i < graph->GetN(); ++i)
				{
					newData[i] = qAdc[i] - Q*exp(-qMeasurementTimes[i]/shapingTime);
				}

				const Double_t newAdcPeak = TMath::MaxElement(numEntries,newData);
				const Double_t newTimePeak = TMath::LocMax(numEntries,newData);

				fittingFunction = new TF1("fittingFunction",fittingFunction7Uniform,firstBin,lastBin,5);
				const double timeShift = newTimePeak - shapingTime;
				const double scalingFactor = TMath::Max(newAdcPeak /0.015, 1000.0);
				const double sigma = 10.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,Q,sigma,shapingTime);
				//fittingFunction->SetParLimits(2,0.0,1000.0);
				//fittingFunction->SetParLimits(0,0.0-shapingTime,140.0-shapingTime); //This is a concession at 0.0
				//fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				//fittingFunction->SetParLimits(3,0.0,30.0);
				fittingFunction->FixParameter(4,shapingTime);
			}
			if (tPeak.size() == 2)
			{
				option = 2;
				fittingFunction = new TF1("fittingFunction",fittingFunction7Uniform,firstBin,lastBin,5);
				const double timeShift = tPeak[1] - shapingTime;
				const double scalingFactor = TMath::Max(adcPeak[1] /0.015, 1000.0);
				const double Q = TMath::Max(qAdc[0],0.0);
				const double sigma = 10.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,Q,sigma,shapingTime);
				//fittingFunction->SetParLimits(2,0.0,1000.0);
				//fittingFunction->SetParLimits(0,0.0-shapingTime,140.0-shapingTime); //This is a concession at 0.0
				//fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				//fittingFunction->SetParLimits(3,0.0,30.0);
				fittingFunction->FixParameter(4,shapingTime);
			}
			if (tPeak.size() == 3)
			{
				option = 3;
				fittingFunction = new TF1("fittingFunction", fittingFunction8,firstBin,lastBin, 6);
				const double timeShift0 = tPeak[1] - shapingTime;
				const double scalingFactor0 = TMath::Max(adcPeak[1] /0.015, 1000.0);
				const double Q = TMath::Max(qAdc[0],0.0);
				const double timeShift1 = tPeak[2] - tPeak[1];
				const double scalingFactor1 = TMath::Max(adcPeak[2] /0.015, 1000.0);

				fittingFunction->SetParameters(timeShift0, scalingFactor0, Q, timeShift1, scalingFactor1, shapingTime);
				//fittingFunction->SetParLimits(0,0.0-shapingTime,80.0-shapingTime);
				//fittingFunction->SetParLimits(1,1000.0,1.0e9);
				//fittingFunction->SetParLimits(3,25.0,80.0);
				//fittingFunction->SetParLimits(4,1000.0,1.0e9);
				fittingFunction->FixParameter(5,shapingTime); 
			}
		}
		else
		{
			if (tPeak.size() == 1)
			{
				option = 4;
				fittingFunction = new TF1("fittingFunction",fittingFunction4Uniform,firstBin,lastBin,6);
				const double timeShift = tPeak[0] - shapingTime;
				const double scalingFactor = TMath::Max(adcPeak[0] /0.015, 1000.0);
				const double verticalShift = 0.5 * (qAdc[0] + qAdc[1]);
				const double sigma = 10.0;
				const double truncationLevel = 1024.0 - 64.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
				//fittingFunction->SetParLimits(0,20.0-shapingTime,140.0-shapingTime); 
				//fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				//fittingFunction->SetParLimits(2,-20.0,1000.0); // This bound currently works but is questionable
				//fittingFunction->SetParLimits(3,0.0,30.0);
				//fittingFunction->FixParameter(4,1024.0-64.0); 
				fittingFunction->FixParameter(5,shapingTime);
			}
			if (tPeak.size() == 2)
			{
				option = 5;
				fittingFunction = new TF1("fittingFunction",fittingFunction10,firstBin,lastBin,6);
				const double timeShift0 = tPeak[0] - shapingTime;
				const double scalingFactor0 = TMath::Max(adcPeak[0] /0.015, 1000.0);
				const double verticalShift = 0.5 * (qAdc[0] + qAdc[1]);
				const double timeShift1 = tPeak[1] - tPeak[0];
				const double scalingFactor1 = TMath::Max(adcPeak[1] /0.015, 1000.0);
				fittingFunction->SetParameters(timeShift0, scalingFactor0, verticalShift, timeShift1, scalingFactor1, shapingTime);
				//fittingFunction->SetParLimits(0,20.0-shapingTime,100.0-shapingTime);
				//fittingFunction->SetParLimits(1,1000.0,1.0e9);
				//fittingFunction->SetParLimits(2,-20.0,1000.0); // This bound currently works but is questionable
				//fittingFunction->SetParLimits(3,25.0,100.0);
				//fittingFunction->SetParLimits(4,1000.0,1.0e9);
				//fittingFunction->FixParameter(5,shapingTime); 
			}
		}

		graph->Fit(fittingFunction,"QN");
		graph->Draw("A*");
		fittingFunction->Draw("same");
//}
	}

