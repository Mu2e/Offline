// Constructor of a PDF to extract random times to describe the proton pulse
//
// $Id: ProtonPulseRandPDF.cc,v 1.18 2014/04/14 19:08:46 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/14 19:08:46 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "Mu2eUtilities/inc/ProtonPulseRandPDF.hh"

// cetlib includes
#include "cetlib/exception.h"

// The following defines the proton pulse shape parameters (pdf width,
// pdf step and differential distribution).  Please note that it is
// preferred to assume a delta function distribution for the proton
// pulse during framework jobs, and then to convolute the obtained
// timing distributions with the desired shape after the fact.
//
// Previous read-in tables (as of Jan. 27, 2014) are deprecated.  
//
// The POT pulse contains two components:
//
//   (1) The central POT pulse, which is not weighted by any effects
//       from the AC dipole
//   (2) The AC dipole transmission function which serves to cut off
//       any out-of-time POTs.
//

namespace mu2e{
  
  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                                         const std::string pulseType ) 
    : accPar_      ( &*ConditionsHandle<AcceleratorParams>( "ignored" ) )
    , pulseShape_  ( loadTable<2>( accPar_->potPulse , false ) )
    , acdipole_    ( loadTable<2>( accPar_->acDipole , false ) )
    , extFactor_   ( accPar_->intrinsicExtinction )
    , pulseEnum_   ( pulseType )
    , times_       ( setTimes() )
    , spectrum_    ( setSpectrum() )
    , randSpectrum_( engine, &spectrum_.front(), times_.size() )
  {
  }  
  

  //============================================================================================================
  double ProtonPulseRandPDF::fire() {

    static const double pdfWidth = timeMax_-timeMin_;
    return pdfWidth*(randSpectrum_.fire() - fireOffset_) ; // subtractive constant centers pulse at 0.
    
  }
  
  //============================================================================================================
  std::vector<double> ProtonPulseRandPDF::setTimes() {
    
    // Determine min. pulse bin width (assume all bins uniformly distributed)
    const double potBinWidth = std::abs( pulseShape_(1,0)-pulseShape_(0,0) );
    const double acBinWidth  = std::abs( acdipole_  (1,0)-acdipole_  (0,0) );

    const double pdfBinWidth = std::min( potBinWidth, acBinWidth );

    if ( pulseEnum_ == PotSpectrumEnum::DEFAULT ) { 
      timeMax_    =  accPar_->limitingHalfWidth;
      timeMin_    = -accPar_->limitingHalfWidth;
    }
    if ( pulseEnum_ == PotSpectrumEnum::TOTAL   ) { 
      timeMax_    =  accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;                        
      timeMin_    = -accPar_->limitingHalfWidth;
    }
    if ( pulseEnum_ == PotSpectrumEnum::OOT     ) { 
      timeMax_    =  accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;
      timeMin_    =  accPar_->limitingHalfWidth;              
    }
    if ( pulseEnum_ == PotSpectrumEnum::ALLFLAT ) {
      timeMax_    =  accPar_->deBuncherPeriod - accPar_->limitingHalfWidth;
      timeMin_    = -accPar_->limitingHalfWidth;
    }

    fireOffset_ =  -timeMin_/(timeMax_-timeMin_);
    const std::size_t nbins = (timeMax_ - timeMin_)/pdfBinWidth+1;

    // Set time axis
    std::vector<double> times;
    times.reserve( nbins );

    for ( std::size_t i(0); i < nbins ; i++ ) 
      times.push_back( timeMin_+i*pdfBinWidth );

    return times;

  };

  //============================================================================================================
  std::vector<double> ProtonPulseRandPDF::setSpectrum() const {

    std::vector<double> spectrum;
    spectrum.reserve( times_.size() );

    if ( pulseEnum_ == PotSpectrumEnum::ALLFLAT ) spectrum = std::vector<double>( times_.size() , 1. );
    else                                          spectrum = preparePotSpectrum();

    return spectrum;
  }

  //============================================================================================================
  std::vector<double> ProtonPulseRandPDF::preparePotSpectrum() const {

    auto potShape    = getShape( pulseShape_ );

    // For convenience, normalize potShape so that portion within limiting halfwidth
    // wrt pulse center has area of 1.
    const double ext = extFactor_; // Once I have full pulse, I can calculate it using: determineIntrinsicExt( potShape, accPar_->limitingHalfWidth );
    renormalizeShape( potShape, 1.+ext );

    // Now replace spectrum outside of halfwidth with flat
    // distribution given an expected intrinsic extinction level
    replaceOotShape( potShape, accPar_->limitingHalfWidth, ext );

    // Get AC dipole transmission function 
    auto dipoleShape  = getShape( acdipole_ ); 

    // Check if maps above have equal key ranges
    if ( potShape.size() != dipoleShape.size() ) 
      throw cet::exception("POTshape") <<
        " POT shape and AC dipole transmission function do not have same set of times!\n ";
    
    std::equal ( potShape.begin(), potShape.end(), dipoleShape.begin(),
                 [](std::pair<const double,double>& pair1,
                    std::pair<const double,double>& pair2){
                   return pair1.first == pair2.first;
                 } );

    // Set the spectrum given shapes above
    std::vector<double> spectrum;
    spectrum.reserve( times_.size() );

    for ( const auto& t : times_ ) 
      spectrum.push_back( potShape.find(t)->second*dipoleShape.find(t)->second );

    return spectrum;

  }

  //============================================================================================================
  std::map<double,double> ProtonPulseRandPDF::getShape( const Table<2>& table, const double timeOffset ) const {

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

    // Linearly interpolate to fill out shape
    std::map<double,double> shape;

    for ( const auto& t : times_ ) {
      if      ( t <= begin->first ) shape[t] = 0.;
      else if ( t >  end->first   ) shape[t] = 0.;
      else {
        auto const& it1 = rawShape.lower_bound( t );
        auto it0 = it1; --it0;
        
        shape[t] = it0->second + (it1->second - it0->second)/(it1->first - it0->first)*(t-it0->first);
      }
    }

    return shape;

  }

  //============================================================================================================
  void ProtonPulseRandPDF::renormalizeShape( std::map<double,double>& shape, const double norm ) {
    
    static auto add_to_integral = [](double integral, const std::pair<const double,double>& it) {
      return integral + it.second;
    };

    const double integral = std::accumulate( shape.begin(), shape.end(), 0., add_to_integral );
    std::for_each( shape.begin(), shape.end(), [&](std::pair<const double,double>& pt){ pt.second *= norm/integral; } );

  }

  //============================================================================================================
  double ProtonPulseRandPDF::determineIntrinsicExt( const std::map<double,double>& shape, const double hw ) {

    double num(0.), denom(0.);
    std::for_each( shape.begin(), shape.end(), 
                   [&](const std::pair<const double,double>& pt ) 
                   {
                     if ( std::abs(pt.first) >= hw ) num += pt.second;
                     else denom += pt.second;
                   } );
    
    return num/denom;
  }

  //============================================================================================================
  void ProtonPulseRandPDF::replaceOotShape( std::map<double,double>& shape, const double hw, const double ext ) {

    unsigned bins(0);
    std::for_each( shape.begin(), shape.end(), [&](const std::pair<const double,double>& pt) { if (std::abs(pt.first) >= hw ) bins++;               } );
    std::for_each( shape.begin(), shape.end(), [&](      std::pair<const double,double>& pt) { if (std::abs(pt.first) >= hw ) pt.second = ext/bins; } );
    
  }

}
