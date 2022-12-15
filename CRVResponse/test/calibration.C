void calibration(std::string filename, int nPEpeaksToFit=1)
{
  TFile *file = new TFile(filename.c_str());
  if(!file->cd("CrvCalibration"))
  {
    std::cout<<"CrvCalibration directory not found in Root file."<<std::endl;
    return;
  }

  const std::string baseName="crvCalibrationHist_";
  TF1 funcCalibPeaks("f1", "gaus");
  TF1 funcCalib("f2","[0]*x");
  for(const auto&& key: *gDirectory->GetListOfKeys())
  {
    TH1 *hist = dynamic_cast<TH1F*>(gDirectory->Get(key->GetName()));
    if(hist==NULL) continue;  //doesn't seem to be a histogram
    const std::string histName=hist->GetName();
    if(histName.compare(0,baseName.length(),baseName)!=0) continue;  //doesn't seem to be a pedestal histogram
    if(hist->GetEntries()<200) continue;  //not enough entries

    int maxbin = 0;
    double maxbinContent = 0;
    int startBin=hist->FindBin(250); //start at 250
    for(int bin=startBin; bin<hist->GetNbinsX(); bin++)
    {
      double binContent = hist->GetBinContent(bin);
      if(binContent>maxbinContent)
      {
        maxbin=bin;
        maxbinContent=binContent;
      }
    }

    //find 1PE peak
    double peak1PE = hist->GetBinCenter(maxbin);
    funcCalibPeaks.SetRange(peak1PE*0.8,peak1PE*1.2);
    funcCalibPeaks.SetParameter(1,peak1PE);
    hist->Fit(&funcCalibPeaks, "NQR");
    peak1PE = funcCalibPeaks.GetParameter(1);

    //find other PE peaks
    std::vector<double> peaks;
    if(peak1PE>250.0 && peak1PE<750.0)
    {
      peaks.push_back(peak1PE);
      for(int iPeak=2; iPeak<=nPEpeaksToFit; ++iPeak)
      {
        double peakTmp = iPeak*peak1PE;
        double rangeStart = peakTmp-0.35*peak1PE;
        double rangeEnd = peakTmp+0.35*peak1PE;
        if(rangeEnd>hist->GetXaxis()->GetXmax()) break;
        if(hist->Integral(hist->FindBin(rangeStart),hist->FindBin(rangeEnd))<50) break;
        funcCalibPeaks.SetRange(rangeStart,rangeEnd);
        funcCalibPeaks.SetParameter(1,peakTmp);
        if((int)hist->Fit(&funcCalibPeaks, "NQR")!=0) break;
        peakTmp = funcCalibPeaks.GetParameter(1);
        if(peakTmp/peak1PE>iPeak*0.95 && peakTmp/peak1PE<iPeak*1.05)
        {
          peaks.push_back(peakTmp);
        }
        else break;
      }
    }
    else continue;

    //plot to determine calibration constants
    TGraph *graph=new TGraph();
    graph->SetPoint(0,0,0);
    for(size_t iPeak=0; iPeak<peaks.size(); ++iPeak) graph->SetPoint(iPeak+1,iPeak+1,peaks[iPeak]);
    funcCalib.SetRange(-0.5, peaks.size()+0.5);
    graph->Fit(&funcCalib, "NQR");
    std::cout<<histName.substr(baseName.length())<<"  "<<funcCalib.GetParameter(0)<<std::endl;
  }
}
