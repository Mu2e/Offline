#include "MakeCrvPhotonArrivals.hh"

#include <TCanvas.h>
#include <TH2F.h>
#include <sstream>
#include <vector>
#include <iostream>

void MakeCrvPhotonArrivals::DrawHistograms()
{
  TCanvas c1("ArrivalProbabilities1","");
  TH2F h1("HistArrivalProbabilities1","",_LBD.yBins.size()-1,_LBD.yBins.data(),_LBD.zBins.size()-1,_LBD.zBins.data());

  for(unsigned int iy=1; iy<_LBD.yBins.size(); iy++)
  for(unsigned int iz=1; iz<_LBD.zBins.size(); iz++)
  {
    double y=(_LBD.yBins[iy-1]+_LBD.yBins[iy])/2.0;
    double z=(_LBD.zBins[iz-1]+_LBD.zBins[iz])/2.0;
    int i=_LBD.findScintillatorBin(0.0,y,z);
    if(i<0) continue;
    const LookupBin &bin = _bins[0][i];
    float p = bin.arrivalProbability[0];
    if(!isnan(p)) h1.Fill(y,z,p);
  }

  h1.Draw("COLZ");
  c1.SaveAs("ArrivalProbability1.C");

  TCanvas c2("ArrivalProbabilities2","");
  std::vector<TH1F*> h2;
  for(double x=-8; x<10; x+=2)
  {
    std::stringstream s;
    s<<"HistArrivalProbabilities2_"<<x+10;
    TH1F *h2Tmp=new TH1F(s.str().c_str(),"",_LBD.zBins.size()-1,_LBD.zBins.data());
    h2Tmp->SetMinimum(0);
    h2Tmp->SetLineWidth(4);

    for(unsigned int iz=1; iz<_LBD.zBins.size(); iz++)
    {
      double z=(_LBD.zBins[iz-1]+_LBD.zBins[iz])/2.0;
      int i=_LBD.findScintillatorBin(x,0.0,z);
      if(i<0) continue;
      const LookupBin &bin = _bins[0][i];
      float p = bin.arrivalProbability[0];
      if(!isnan(p)) h2Tmp->Fill(z,p);
    }
    h2Tmp->Draw("same");
    h2.push_back(h2Tmp);
  }

  c2.SaveAs("ArrivalProbability2.C");
}

