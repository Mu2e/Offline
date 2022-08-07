void digiGain()
{
  std::string filename = "stmDigisSpectrum.root";
  TFile* data = new TFile(filename.c_str(), "READ");
  /*
  // Get subset of adcSpectrum
  double threshold = 500.0;
  std::vector<double> adcSubset;
  for(int i = 0; i < data->Get("plotSTMDigisSpectrum/adcSpectrum").size(); ++i) {
    if (data[i] >= threshold) {
      adcSubset.push_back(data[i]);
    }
  }
  TH1D* adcSpectrum = (TH1D*) adcSubset;
  */
  // Without subset
  TH1D* adcSpectrum = (TH1D*) data->Get("plotSTMDigisSpectrum/adcSpectrum");
  adcSpectrum->Rebin(10);

  // Initializing the energy peak and adc peak locations
  const double energyPeaks[4] = {0.898,1.173,1.333,1.836};
  std::vector<double> adcPeaks;

  // Finding the location of each peak
  const int n_peaks = 8;
  TSpectrum* spectrum = new TSpectrum(n_peaks);
  int n_found_peaks = spectrum->Search(adcSpectrum);
  for (int i_peak = 0; i_peak < n_found_peaks; ++i_peak)
    {
      double peak_x_pos = *(spectrum->GetPositionX() + i_peak);
      std::cout << peak_x_pos << std::endl;
      adcPeaks.push_back(peak_x_pos);
    }
  // Initializing the output file
  ofstream out;
  out.open("adcPeaks.log", ios::out | ios::trunc);
  out << "peak_number, par_name, par_value, par_error" << std::endl;

  // Plotting the adc spectrum
  adcSpectrum->Draw("HIST");

  // Estimating background using TSpectrum
  TH1* hb = spectrum->Background(adcSpectrum, 20, "SAME");

  TF1* fline = new TF1("fline", "pol1", 0, 5000);
  adcSpectrum->Fit(fline, "qn");

  //adcSpectrum->Add(hb,-1);

  //Fitting
  for (int j = 0; j < adcPeaks.size(); ++j)
    {
      TString fname(Form("fgaus_%d", j));
      TF1* fitGaus = new TF1(fname, "[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",adcPeaks[j]-50,adcPeaks[j]+50);
      fitGaus->SetParName(0,"Amplitude");
      fitGaus->SetParName(1,"Mean");
      fitGaus->SetParName(2,"Sigma");
      fitGaus->SetParName(3,"Linear");
      fitGaus->SetParName(4,"Constant");

      // Guessing the height for each TSpectrum peak
      int bin = adcSpectrum->FindBin(adcPeaks[j]);
      double guessHeight = adcSpectrum->GetBinContent(bin);
      std::cout << "TSpectrum value for peak " << j  <<  ": " << int (adcPeaks[j]) << std::endl;
      std::cout << "Guess for height of peak " << j  <<  ": " << guessHeight       << std::endl;

      // Guess for FWHM
      double bkrndHeight = hb->GetBinContent(bin);
      double sgnlHeight = guessHeight - bkrndHeight;
      double halfMax = sgnlHeight/2.0;
      double x1;
      double x2;
      bool threshold = false;
      std::cout << "Half Max: " << halfMax << std::endl;
      for (int k = bin - 5; k < bin + 5; ++k)
        {
          double k_content = adcSpectrum->GetBinContent(k) - hb->GetBinContent(k);
          //std::cout << "k_content: " << k_content << std::endl;
          if (k_content > halfMax && threshold == false)
            {
              x1 = adcSpectrum->GetBinCenter(k);
              threshold = true;
              //std::cout << "x1: " << x1 << std::endl;
            }
          if (k_content < halfMax && threshold == true)
            {
              x2 = adcSpectrum->GetBinCenter(k);
              //std::cout << "x2: " << x2 << std::endl;
              break;
            }
        }
      double FWHM = x2 - x1;
      //std::cout << "FWHM: " << FWHM << std::endl;
      double sigma = FWHM/2.35;
      std::cout << "Sigma: " << sigma << std::endl;

      fitGaus->SetParameters(guessHeight,adcPeaks[j],sigma,-0.05,2);
      fitGaus->SetParLimits(1,adcPeaks[j]-10,adcPeaks[j]+10);
      //fitGaus->SetParLimits(0,0,200);
      //fitGaus->SetParLimits(2,0.001,10);

      // If the found peak is less than 500 ADC sample then don't fit
      if(adcPeaks[j] > 500) {
      TFitResultPtr fitresult = adcSpectrum->Fit(fitGaus, "RS+");
      int n_par = fitresult->NPar();
      //Draw fitted peaks
      fitGaus->Draw("LSAME");
      // Export fit parameters to a file
      for (int i_par = 0; i_par < n_par; ++i_par)
        {
          out   << j << ", "
                << fitresult->GetParameterName(i_par)
                << ", " << fitresult->Parameter(i_par)
                << ", " << fitresult->ParError(i_par)
                << std::endl;
        }
      }
    }
  out.close();
}
