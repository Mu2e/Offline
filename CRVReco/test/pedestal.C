void output(std::ofstream &outputFile, int channel, double pedestal)
{
  outputFile<<channel<<","<<pedestal<<",-1,-1"<<std::endl;
}

void pedestal(const std::string &filename)
{
  TFile *file = new TFile(filename.c_str());
  if(!file->cd("CrvPedestalFinder"))
  {
    std::cout<<"CrvPedestalFinder directory not found in Root file."<<std::endl;
    return;
  }

  std::filesystem::path path(filename);
  path.replace_extension("txt");
  std::ofstream outputFile;
  outputFile.open(path);

  /***************************************************************************/

  //CRVSiPM table - only fill the pedestal, set all other values to 0
  outputFile<<"TABLE CRVSiPM"<<std::endl;
  outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

  const std::string baseName="crvPedestalHist_";
  TF1 funcPedestal("f0", "gaus");
  for(const auto&& key: *gDirectory->GetListOfKeys())
  {
    TH1 *hist = dynamic_cast<TH1F*>(gDirectory->Get(key->GetName()));
    if(hist==NULL) continue;  //doesn't seem to be a histogram
    const std::string histName=hist->GetName();
    if(histName.compare(0,baseName.length(),baseName)!=0) continue;  //doesn't seem to be a pedestal histogram
    const std::string channelName=histName.substr(baseName.length());
    int channel=atoi(channelName.c_str());

    if(hist->GetEntries()<200) {output(outputFile,channel,0); continue;}  //not enough entries
    int n=hist->GetNbinsX();
    double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
    if(overflow/((double)hist->GetEntries())>0.1) {output(outputFile,channel,0); continue;}  //too much overflow. something may be wrong

    int maxbinPedestal = hist->GetMaximumBin();
    double peakPedestal = hist->GetBinCenter(maxbinPedestal);
    funcPedestal.SetRange(peakPedestal-4,peakPedestal+4);
    funcPedestal.SetParameter(1,peakPedestal);
    hist->Fit(&funcPedestal, "NQR");
    output(outputFile,channel,funcPedestal.GetParameter(1));
  }

  outputFile<<std::endl;

  /***************************************************************************/

  //CRVTime table - set all values to 0
  outputFile<<"TABLE CRVTime"<<std::endl;
  outputFile<<"#channel, timeOffset"<<std::endl;
  for(const auto&& key: *gDirectory->GetListOfKeys())
  {
    TH1 *hist = dynamic_cast<TH1F*>(gDirectory->Get(key->GetName()));
    if(hist==NULL) continue;  //doesn't seem to be a histogram
    const std::string histName=hist->GetName();
    if(histName.compare(0,baseName.length(),baseName)!=0) continue;  //doesn't seem to be a pedestal histogram
    const std::string channelName=histName.substr(baseName.length());
    outputFile<<channelName<<",0"<<std::endl;
  }


  /***************************************************************************/

  outputFile.close();
  gApplication->Terminate();
}
