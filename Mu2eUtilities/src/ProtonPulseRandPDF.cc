//
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.7 2013/03/28 17:21:11 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/03/28 17:21:11 $
//
// Original author Gianni Onorato
//

// Framework includes
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

//C++ includes
#include <cmath>
#include <fstream>
//#include <utility>
//CLHEP includes
#include "CLHEP/Random/RandFlat.h"

// The following defines the proton pulse shape parameters (pdf width,
// pdf step and differential distribution).  Please note that it is
// preferred to assume a delta function distribution for the proton
// pulse during framework jobs, and then to convolute the obtained
// timing distributions with the desired shape after the fact.
//
// The shape below supersedes "ProtonPulseSpectrum.txt" and is current
// as of Mar. 2013.  Also note that in any framework jobs that have
// the proton pulse shape enabled, the center of the pulse will be
// shifted 200 ns wrt 0.

using namespace std;

static const double _PDFwidth = 400; //time in nanosecond
static const double _PDFstep  = 2.; 
static const string _FILEname = "ConditionsService/data/ProtonPulseSpectrum_01.txt";

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
    return nbins;
  }
  
  vector<double> ProtonPulseRandPDF::ProtonPulseSpectrum() {

    int count = 0;
    
    vector<double> spectrum;
    ConfigFileLookupPolicy findConfig;
    const string infile = findConfig( _FILEname );
    fstream inst(infile.c_str(), ios::in);
    if (!inst.is_open()) {
      throw cet::exception("ProductNotFound")
	<< "No proton spill input file found";
    }
    double invalue;
    while (!(inst.eof())) {
      inst >> invalue;
      count++;
      spectrum.push_back(invalue);
    }
    return spectrum;
  }
}
