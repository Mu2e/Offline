void pedestal(std::string filename)
{
  TFile *file = new TFile(filename.c_str());
  if(!file->cd("CrvPedestalFinder"))
  {
    std::cout<<"CrvPedestalFinder directory not found in Root file."<<std::endl;
    return;
  }

  const std::string baseName="crvPedestalHist_";
  TF1 funcPedestal("f0", "gaus");
  for(const auto&& key: *gDirectory->GetListOfKeys())
  {
    TH1 *hist = dynamic_cast<TH1F*>(gDirectory->Get(key->GetName()));
    if(hist==NULL) continue;  //doesn't seem to be a histogram
    const std::string histName=hist->GetName();
    if(histName.compare(0,baseName.length(),baseName)!=0) continue;  //doesn't seem to be a pedestal histogram
    if(hist->GetEntries()<200) continue;  //not enough entries

    int n=hist->GetNbinsX();
    double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
    if(overflow/((double)hist->GetEntries())>0.1) continue;  //too much overflow. something may be wrong

    int maxbinPedestal = hist->GetMaximumBin();
    double peakPedestal = hist->GetBinCenter(maxbinPedestal);
    funcPedestal.SetRange(peakPedestal-4,peakPedestal+4);
    funcPedestal.SetParameter(1,peakPedestal);
    hist->Fit(&funcPedestal, "NQR");
    std::cout<<histName.substr(baseName.length())<<"  "<<funcPedestal.GetParameter(1)<<std::endl;
  }
}
