// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.16 2014/04/01 15:03:16 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/01 15:03:16 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

// cetlib includes
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
//   (1) The central POT pulse, which is not weighted by any effects
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
  const double ootOffset = 500. ; // ns 
}

namespace mu2e{
  
  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                                         const std::string pulseType ) 
    : accPar_      ( &*ConditionsHandle<AcceleratorParams>( "ignored" ) )
    , pulseShape_  ( loadTable<2>( accPar_->potPulse , false ) )
    , acdipole_    ( loadTable<2>( accPar_->acDipole , false ) )
    , ootPulse_    ( loadTable<2>( accPar_->outOfTime, false ) )
    , pulseEnum_   ( pulseType )
    , nPoints_     ( calculateNpoints() )
    , extFactor_   ( 1e-10 )
    , spectrum_    ( setSpectrum() )
    , randSpectrum_(engine, &spectrum_.front(), nPoints_ )
  {
  }  
  
  double ProtonPulseRandPDF::fire() {

    static const double pdfWidth = timeMax_-timeMin_;
    return pdfWidth*(randSpectrum_.fire() - fireOffset_) ; // subtractive constant centers pulse at 0.
    
  }
  
  std::size_t ProtonPulseRandPDF::calculateNpoints() {
    
    ConditionsHandle<AcceleratorParams> accPar ("ignored");

    // Determine min. pulse bin width (assume all bins uniformly distributed)
    const double potBinWidth = std::abs( pulseShape_(1,0)-pulseShape_(0,0) );
    const double acBinWidth  = std::abs( acdipole_  (1,0)-acdipole_  (0,0) );
    const double ootBinWidth = std::abs( ootPulse_  (1,0)-ootPulse_  (0,0) );

    const double pdfBinWidth = std::min( std::min( potBinWidth, acBinWidth ), ootBinWidth );

    if ( pulseEnum_ == PotSpectrum::DEFAULT ) { 
      timeMax_    = accPar_->limitingHalfWidth;
      timeMin_    = -accPar_->limitingHalfWidth;
    }
    if ( pulseEnum_ == PotSpectrum::TOTAL   ) { 
      timeMax_    = accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;                        
      timeMin_    = pulseShape_(0,0);
    }
    if ( pulseEnum_ == PotSpectrum::OOT     ) { 
      timeMax_    = accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;
      timeMin_    = accPar_->limitingHalfWidth;              
    }
    if ( pulseEnum_ == PotSpectrum::OOTFLAT ) {
      timeMax_    = accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;
      timeMin_    = accPar_->limitingHalfWidth;
    }
    if ( pulseEnum_ == PotSpectrum::ALLFLAT ) {
      timeMax_    = accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;
      timeMin_    = -accPar_->limitingHalfWidth;
    }

    fireOffset_ =  -timeMin_/(timeMax_-timeMin_);
    const std::size_t nbins = (timeMax_ - timeMin_)/pdfBinWidth+1;

    // Set time axis
    times_.reserve( nbins );
    for ( std::size_t i(0); i < nbins ; i++ ) 
      times_.push_back( timeMin_+i*pdfBinWidth );

    return nbins;

  };


  std::vector<double> ProtonPulseRandPDF::getShape( const Table<2>& table, const double timeOffset ) const {

    // Get raw shape
    std::map<double,double> rawShape;

    for ( std::size_t i(0) ; i < table.getNrows() ; i++ ) {
      const double time = table(i,0) - timeOffset;
      rawShape.insert( std::make_pair( time, table(i,1) ) );  
      if ( time > timeMax_ ) break;
    }

   // Get shape boundaries
    auto const& begin = rawShape.begin();
    auto end = rawShape.end(); --end;

    // Interpolate to fill out shape
    std::vector<double> shape;
    shape.reserve( nPoints_ );

    for ( std::size_t i(0) ; i < nPoints_ ; i++ ) {
      if      ( times_.at(i) <= begin->first ) shape.push_back(0.);
      else if ( times_.at(i) >  end->first   ) shape.push_back(0.);
      else {
        auto const& it1 = rawShape.lower_bound( times_.at(i) );
        auto it0 = it1; --it0;
        
        const double intSpectrum = it0->second + (it1->second - it0->second)/(it1->first - it0->first)*(times_.at(i)-it0->first);
        shape.push_back( intSpectrum );
      }
    }

    return shape;
  }

  std::vector<double> ProtonPulseRandPDF::setSpectrum() const {

    std::vector<double> potShape    = getShape( pulseShape_ );

    // For convenience, normalize potShape to unity
    renormalizeShape( potShape, 1. );

    // Get AC dipole transmission function and out-of-time shape
    std::vector<double> dipoleShape  = getShape( acdipole_   ); 
    std::vector<double> ootShape     = getShape( ootPulse_, ootOffset );
    std::vector<double> flatShape ( nPoints_, 1. ); 

    // Normalize out-of-time shapes given extinction factor
    renormalizeShape( ootShape  , extFactor_ );
    renormalizeShape( flatShape , extFactor_ );

    // Set the spectrum given shapes above
    std::vector<double> spectrum;
    spectrum.reserve( nPoints_ );

    for ( std::size_t i(0) ; i < nPoints_ ; i++ ) { 
      double weight(0.);
      if ( pulseEnum_ == PotSpectrum::TOTAL   ) weight = potShape.at(i)*dipoleShape.at(i) + ootShape.at(i);
      if ( pulseEnum_ == PotSpectrum::DEFAULT ) weight = potShape.at(i)*dipoleShape.at(i);
      if ( pulseEnum_ == PotSpectrum::OOT     ) weight = ootShape.at(i);
      if ( pulseEnum_ == PotSpectrum::OOTFLAT ) weight = flatShape.at(i);
      if ( pulseEnum_ == PotSpectrum::ALLFLAT ) weight = flatShape.at(i);
      spectrum.push_back( weight );
    }
    
    return spectrum;

  }


  void ProtonPulseRandPDF::renormalizeShape( std::vector<double>& shape, const double norm ) const {
    
    const double integral = std::accumulate( shape.begin(), shape.end(), 0. );
    std::for_each( shape.begin(), shape.end(), [&](double& pt){ pt *= norm/integral; } );

  }
  
}


