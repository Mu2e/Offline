//
// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.1 2011/06/13 17:07:04 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/06/13 17:07:04 $
//
// Original author Gianni Onorato
//

// Framework includes
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

//C++ includes
#include <cmath>

//CLHEP includes
#include "CLHEP/Random/RandFlat.h"


namespace mu2e{

  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine):
    _randFlat(engine)
  {
    binWidth = 0.25; //Precision of the PDF in nanoseconds
    pulseRange = 260.0000; //Width of the proton pulse time distribution  
    sigma1 = 26.5000;
    sigma2 = 25.5000;
    sigma3 = 6.4000;
    A1 = 426.0000; 
    A2 = 433.0000;
    A3 = 252.0000;
    x01 = 130.0000-36.0000;
    x02 = 130.0000+38.0000;
    x03 = 130.0000+0.5300;
    Norm = 60007.17140936479;
    // Calculated with Mathematica from 0 to 260. Above the complete spectrum should be sqrt(CLHEP::twopi)*(A1*sigma1+A2*sigma2+A3*sigma3);
  }
  
  ProtonPulseRandPDF::~ProtonPulseRandPDF()
  {}
  
  double ProtonPulseRandPDF::TripleGaussian(double x) {
    
    double ex1 = (x-x01)*(x-x01)/(2*sigma1*sigma1);
    double G1 = A1 * (exp(-ex1));
    double ex2 = (x-x02)*(x-x02)/(2*sigma2*sigma2);
    double G2 = A2 * (exp(-ex2));
    double ex3 = (x-x03)*(x-x03)/(2*sigma3*sigma3);
    double G3 = A3 * (exp(-ex3));


    return (G1+G2+G3)/Norm;
  }

  double ProtonPulseRandPDF::fire() {

    double Yflat = _randFlat.fire();
    double partialint = 0;
    double toextract = 0;
    while (partialint < Yflat && toextract <= pulseRange) {
      partialint += TripleGaussian(toextract)*binWidth;
      toextract += binWidth;
    }
    return toextract;
  }

}
