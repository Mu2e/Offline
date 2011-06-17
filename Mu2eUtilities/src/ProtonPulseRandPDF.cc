//
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.4 2011/06/17 21:08:27 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/06/17 21:08:27 $
//
// Original author Gianni Onorato
//

// Framework includes
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"

//C++ includes
#include <cmath>
#include <fstream>
//#include <utility>
//CLHEP includes
#include "CLHEP/Random/RandFlat.h"

using namespace std;

static const double _PDFwidth = 260; //time in nanosecond 
static const double _PDFstep  = 0.5; 

namespace mu2e{
  
  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine):
    _nBinsSpectrum( calculateNBins() ),
    _randSpectrum(engine, &(ProtonPulseSpectrum()[0]), _nBinsSpectrum )
  {
  }   
  
  ProtonPulseRandPDF::~ProtonPulseRandPDF()
  {}
  
  double ProtonPulseRandPDF::fire() {
    
    return _PDFwidth * _randSpectrum.fire();
    
  }
  
  int ProtonPulseRandPDF::calculateNBins() {
    
    int nbins = int(_PDFwidth / _PDFstep);
    cout << "NUMBER OF BINS ***********" << nbins << endl;
    return nbins;
  }
  
  vector<double> ProtonPulseRandPDF::ProtonPulseSpectrum() {

    int count = 0;
    
    vector<double> spectrum;
    ConfigFileLookupPolicy findConfig;
    const string infile = findConfig("ConditionsService/data/ProtonPulseSpectrum.txt");
    fstream inst(infile.c_str(), ios::in);
    if (!inst.is_open()) {
      throw cet::exception("ProductNotFound")
	<< "No proton spill input file found";
    }
    double invalue;
    while (!(inst.eof())) {
      inst >> invalue;
      cout << invalue << '\t' << count << endl;
      count++;
      spectrum.push_back(invalue);
    }
    return spectrum;
  }
}
