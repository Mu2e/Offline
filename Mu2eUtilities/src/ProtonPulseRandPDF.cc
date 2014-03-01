// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.13 2014/03/01 18:56:01 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/03/01 18:56:01 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

// Framework includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "cetlib/exception.h"

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
// The POT pulse contains three components:
//
//   (1) The central POT pulse, which is notweighted by any effects
//       from the AC dipole
//   (2) The AC dipole transmission function which cuts off the POT
//       pulse outside of +/- 130 ns
//   (3) The out-of-time distribution of protons that are not removed
//       by the AC dipole.  These contributions are normalized to a given
//       extinction factor.  The current default for the factor is 1e-10,
//       which is explicitly hard-coded into the intialization list below.
//
// Various constants are currently included in the anonymous namespace
// directly below.  Eventually, these need to be included in the
// conditions database.  Note that the out-of-time table is currently
// specified relative to the time the proton reaches the AC dipole.
// The distance between the dipole and the tungsten target corresponds
// to a timing offset of roughly 500 ns.

namespace {

  const double ootCycle  = 1700.; // ns
  const double ootOffset = 500. ; // ns
  const double potWidth  = 260. ; // ns

}

namespace mu2e{
  
  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                                         const std::string pulseType ):
    _pulseShape ( loadTable<2>( "ConditionsService/data/ProtonPulseSpectrum_02.txt"      , false ) ),
    _acdipole   ( loadTable<2>( "ConditionsService/data/ACdipoleTransmissionFunction.txt", false ) ),
    _ootPulse   ( loadTable<2>( "ConditionsService/data/OutOfTimeSpectrum.txt"           , false ) ),
    _pulseEnum  ( pulseType ),
    _nPoints    ( calculateNpoints() ),
    _extFactor  ( 1e-10 ),
    _spectrum   ( setSpectrum() ),
    _randSpectrum(engine, &_spectrum.front(), _nPoints )
  {
  }  
  
  double ProtonPulseRandPDF::fire() {

    static const double pdfWidth = _timeMax-_timeMin;
    return pdfWidth*(_randSpectrum.fire() - _fireOffset) ; // subtractive constant centers pulse at 0.
    
  }
  
  std::size_t ProtonPulseRandPDF::calculateNpoints() {
    
    // Determine min. pulse bin width (assume all bins uniformly distributed)
    const double potBinWidth = std::abs( _pulseShape(1,0)-_pulseShape(0,0) );
    const double acBinWidth  = std::abs( _acdipole  (1,0)-_acdipole  (0,0) );
    const double ootBinWidth = std::abs( _ootPulse  (1,0)-_ootPulse  (0,0) );

    const double pdfBinWidth = std::min( std::min( potBinWidth, acBinWidth ), ootBinWidth );

    if ( _pulseEnum == PotSpectrum::DEFAULT ) { 
      _timeMax    = _pulseShape( _pulseShape.getNrows()-1, 0 ); 
      _timeMin    = _pulseShape(0,0);
      _fireOffset = 0.5; // Pulse should start at beginning of POT pulse
    }
    if ( _pulseEnum == PotSpectrum::TOTAL   ) { 
      _timeMax    = ootCycle-potWidth/2;                        
      _timeMin    = _pulseShape(0,0);
      _fireOffset = -_timeMin/(_timeMax-_timeMin); // Pulse should start at beginning of POT pulse
    }
    if ( _pulseEnum == PotSpectrum::OOT     ) { 
      _timeMax    = ootCycle-potWidth/2;                        
      _timeMin    = 0.;              
      _fireOffset = 0.; // Pulse should start at t = 0. 
    }

    const std::size_t nbins = (_timeMax - _timeMin)/pdfBinWidth+1;

    // Set time axis
    _times.reserve( nbins );
    for ( std::size_t i(0); i < nbins ; i++ ) 
      _times.push_back( _timeMin+i*pdfBinWidth );

    return nbins;

  };


  std::vector<double> ProtonPulseRandPDF::getShape( const Table<2>& table, const double timeOffset ) const {

    // Get raw shape
    std::map<double,double> rawShape;

    for ( std::size_t i(0) ; i < table.getNrows() ; i++ ) {
      const double time = table(i,0) - timeOffset;
      rawShape.insert( std::make_pair( time, table(i,1) ) );  
      if ( time > _timeMax ) break;
    }

   // Get shape boundaries
    auto const& begin = rawShape.begin();
    auto end = rawShape.end(); --end;

    // Interpolate to fill out shape
    std::vector<double> shape;
    shape.reserve( _nPoints );

    for ( std::size_t i(0) ; i < _nPoints ; i++ ) {
      if      ( _times.at(i) <= begin->first ) shape.push_back(0.);
      else if ( _times.at(i) >  end->first   ) shape.push_back(0.);
      else {
        auto const& it1 = rawShape.lower_bound( _times.at(i) );
        auto it0 = it1; --it0;
        
        const double intSpectrum = it0->second + (it1->second - it0->second)/(it1->first - it0->first)*(_times.at(i)-it0->first);
        shape.push_back( intSpectrum );
      }
    }

    return shape;
  }

  std::vector<double> ProtonPulseRandPDF::setSpectrum() const {

    std::vector<double> potShape    = getShape( _pulseShape );

    // For convenience, normalize potShape to unity
    renormalizeShape( potShape, 1. );

    // Get AC dipole transmission function and out-of-time shape
    std::vector<double> dipoleShape = getShape( _acdipole   ); 
    std::vector<double> ootShape    = getShape( _ootPulse, ootOffset );

    // Normalize ootShape given extinction factor
    renormalizeShape( ootShape, _extFactor );

    // Set the spectrum given shapes above
    std::vector<double> spectrum;
    spectrum.reserve( _nPoints );

    for ( std::size_t i(0) ; i < _nPoints ; i++ ) { 
      double weight(0.);
      if ( _pulseEnum == PotSpectrum::TOTAL   ) weight = potShape.at(i)*dipoleShape.at(i) + ootShape.at(i);
      if ( _pulseEnum == PotSpectrum::DEFAULT ) weight = potShape.at(i)*dipoleShape.at(i);
      if ( _pulseEnum == PotSpectrum::OOT     ) weight = ootShape.at(i);
      spectrum.push_back( weight );
    }
    
    return spectrum;

  }


  void ProtonPulseRandPDF::renormalizeShape( std::vector<double>& shape, const double norm ) const {
    
    double integral(0.);
    std::for_each( shape.begin(), shape.end(), [&](double& pt){ integral += pt;      } );
    std::for_each( shape.begin(), shape.end(), [&](double& pt){ pt *= norm/integral; } );

  }
  
}


