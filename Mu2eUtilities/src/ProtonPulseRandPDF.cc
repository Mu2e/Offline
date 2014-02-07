// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.9 2014/02/07 14:48:44 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/07 14:48:44 $
//
// Original author: Gianni Onorato
//                  Kyle Knoepfel (significant updates)

// Mu2e includes
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

// Framework includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

// C++ includes
#include <iostream>
#include <cmath>

// The following defines the proton pulse shape parameters (pdf width,
// pdf step and differential distribution).  Please note that it is
// preferred to assume a delta function distribution for the proton
// pulse during framework jobs, and then to convolute the obtained
// timing distributions with the desired shape after the fact.
//
// Previous read-in tables (as of Jan. 27, 2014) are deprecated, as
// the pdf width and stepsize needed to be explicitly stated because
// the tables were only one column wide.  
//
// Please note that no pdf interpolation is currently supported for
// the proton time pulse.

namespace mu2e{
  
  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine):
    _pulseShape ( loadTable<2>( "ConditionsService/data/ProtonPulseSpectrum_02.txt", false ) ),
    _acdipole   ( loadTable<2>( "ConditionsService/data/ACdipoleTransmissionFunction.txt"  ) ),
    _spectrum   ( setSpectrum() ),
    _timeMin    ( _pulseShape(0,0) ),
    _timeMax    ( _pulseShape( _pulseShape.getNrows()-1, 0 ) ),
    _randSpectrum(engine, &_spectrum.front(), _pulseShape.getNrows() )
  {
  }  
  
  double ProtonPulseRandPDF::fire() {

    static const double pdfWidth = _timeMax-_timeMin;
    return pdfWidth*(_randSpectrum.fire() - 0.5) ; // subtractive constant centers pulse at 0.
    
  }
  
  std::vector<double> ProtonPulseRandPDF::setSpectrum() {

    std::vector<double> shape;
    for ( unsigned iRow(0) ; iRow < _pulseShape.getNrows() ; iRow++ ) { 
      const double time          = _pulseShape( iRow, 0 );
      const unsigned acdipoleRow = std::abs( time ) > 200. ? 0 : _acdipole.findLowerBoundRow( time );

      const double pulseWeight    = _pulseShape( iRow        );
      const double acdipoleWeight = _acdipole  ( acdipoleRow );
 
      shape.push_back( pulseWeight*acdipoleWeight );
    }
    
    return shape;

  }
}
