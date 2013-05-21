#include <vector>
#include <fstream>
#include <iostream>

#include "TH1.h"

#include "ConvoluteHistogram.h"

TH1* ConvoluteHistogram::GetConvHistogram( const TH1& histo, const TH1& hConv ) {
  TH1* convHisto = dynamic_cast<TH1*>(histo.Clone());
  convHisto->Reset();
  
  for ( int i(1) ; i <= histo.GetNbinsX() ; i++ ) {
    const double mean ( histo.GetBinCenter (i) );
    const double area ( histo.GetBinContent(i) );
    for ( int j(0) ; j < _nSamplings ; j++ ) {
      convHisto->Fill( mean + hConv.GetRandom(), area/_nSamplings );
    }
  }
  
  return convHisto;
  
}

TH1* ConvoluteHistogram::GetConvHistogram( const TH1& histo, 
                                           const std::string functionLoc,
                                           const double halfwidth ) {
  std::vector<double> pulse;

  ifstream inf( functionLoc );
  double number;
  int    nentries(0);
  while ( inf >> number ) { pulse.push_back( number ) ; nentries++; }

  TH1D hConv("convolutionShape","Convolution Shape", nentries, -halfwidth, halfwidth );
  
  int iCount(1);
  for ( double const& value : pulse ) {
    hConv.SetBinContent(iCount,value); iCount++;
  }

  return GetConvHistogram( histo, hConv );

}



