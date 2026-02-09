// Constructor of a PDF to extract random times to describe the proton pulse
//
//
// Original author: Kyle Knoepfel

#include <exception>
#include <stdlib.h>
// C++ includes
#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

#include "CLHEP/Random/RandGeneral.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "Offline/Mu2eUtilities/inc/ProtonPulseRandPDF.hh"
#include "Offline/Mu2eUtilities/inc/Table.hh"

// The following defines the proton pulse shape parameters (pdf width,
// pdf step and differential distribution).  Please note that it is
// preferred to assume a delta function distribution for the proton
// pulse during framework jobs, and then to convolute the obtained
// timing distributions with the desired shape directly before
// hit-making is performed.
//
// Read-in tables used before Jan. 27, 2014 are incompatible with the
// current version of this class.
//
// The POT pulse contains two components:
//
//   (1) The central POT pulse, which is not weighted by any effects
//       from the AC dipole
//   (2) The AC dipole transmission function which serves to cut off
//       any out-of-time POTs.
//
// The pulseEnum_ value specifies a shape and a corresponding range of
// sample times, as defined in the setTmin() and setTmax() private
// member functions.  This range can be overridden by explicitly
// specifying values for "tmin" and "tmax" in the FHiCL ParameterSet
// object that is passed as an argument to this class.  One can also
// specify a desired resolution for time range.
//
// A linear interpolation algorithm in the Table<N> class is used if a
// time is specified that does not correspond to one of the loaded
// table entries.
//
// =NB1= No wrapping functionality is available if a time range is
//       specified that is larger than the nominal 1695-ns microbunch
//       time cycle.
//
// =NB2= The proton pulse is always centered at 0 ns.
//
// =NB3= As currently structured, there is an erroneous 1-ns offset in
//       the POT pulse as taken from the potPulse file.  This is
//       because the keys specified in the file correspond to the low-edge
//       value of histogram bins.  This should be fixed in future versions,
//       but is unlikely to result in any noticeable effect.

namespace mu2e{

  ProtonPulseRandPDF::ProtonPulseRandPDF(art::RandomNumberGenerator::base_engine_t& engine,
                                         const Config& conf)
    : pulseEnum_   ( conf.pulseType() )
    , limitingHalfWidth_  ( conf.limitingHalfWidth() )
    , DRPeriod_    ( GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod() )
    , tmin_        ( setTmin() )
    , tmax_        ( setTmax() )
    , tres_        ( conf.tres() )
    , times_       ( setTimes() )
    , extFactor_   ( determineIntrinsicExt( ConfigFileLookupPolicy()( conf.potPulse() ) ) )
    , pulseShape_  ( setPotPulseShape(      ConfigFileLookupPolicy()( conf.potPulse() ) ) )
    , acdipole_    ( loadTable<2>(          ConfigFileLookupPolicy()( conf.acDipole() ) ) )
    , spectrum_    ( setSpectrum() )
    , randSpectrum_( engine, &spectrum_.front(), times_.size() )
  {
    // Process optional FHICL parameters: modify the values set in the
    // initialization list if job config supplies overrides.
    conf.tmin(tmin_);
    conf.tmax(tmax_);
  }


  //============================================================================================================
  double ProtonPulseRandPDF::fire() {
    static const double pdfWidth   = tmax_-tmin_;
    static const double fireOffset = -tmin_/(tmax_-tmin_);
    return pdfWidth*(randSpectrum_.fire() - fireOffset) ; // subtractive constant centers pulse at 0.
  }

  //============================================================================================================
  double ProtonPulseRandPDF::setTmin() {
    double min(0.);
    if ( pulseEnum_ == DEFAULT ) min = -limitingHalfWidth_;
    if ( pulseEnum_ == TOTAL   ) min = -limitingHalfWidth_;
    if ( pulseEnum_ == OOT     ) min =  limitingHalfWidth_;
    if ( pulseEnum_ == ALLFLAT ) min = -limitingHalfWidth_;
    return min;
  }

  //============================================================================================================
  double ProtonPulseRandPDF::setTmax() {
    double max(0.);
    if ( pulseEnum_ == DEFAULT ) max = limitingHalfWidth_;
    if ( pulseEnum_ == TOTAL   ) max = DRPeriod_ - limitingHalfWidth_;
    if ( pulseEnum_ == OOT     ) max = DRPeriod_ - limitingHalfWidth_;
    if ( pulseEnum_ == ALLFLAT ) max = DRPeriod_ - limitingHalfWidth_;
    return max;
  }

  //============================================================================================================
  std::vector<double> ProtonPulseRandPDF::setTimes() {
    std::vector<double> times;
    for ( double t = tmin_ ; t <= tmax_ ; t += tres_ ) times.push_back( t );
    return times;
  }

  //============================================================================================================
  TableVec<2> ProtonPulseRandPDF::setPotPulseShape( const std::string& shapeTxtFile ) {

    Table<2> pulseShape = loadTable<2>( shapeTxtFile );

    // For convenience, normalize potShape so that portion within limiting halfwidth
    // wrt pulse center has area of 1.
    pulseShape.renormalizeShape( 1.+extFactor_ );

    // Create new vector of times, using interpolated values
    const double potBinWidth  = std::abs( pulseShape(1,0)-pulseShape(0,0) );
    const double timeBinWidth = std::abs( times_.at(1)-times_.at(0) );
    const double binCorrectionFactor = timeBinWidth/potBinWidth;

    TableVec<2> newshapevec;
    for ( const double& t : times_ )
      newshapevec.emplace_back( t, pulseShape.getValueAtKey(t, binCorrectionFactor ) );

    // Now replace spectrum outside of halfwidth with flat
    // distribution given an expected intrinsic extinction level
    unsigned bins(0);
    std::for_each( newshapevec.begin(), newshapevec.end(), [&](const TableRow<2>& pt) { if (std::abs(pt.first) >= limitingHalfWidth_ ) ++bins;                            } );
    std::for_each( newshapevec.begin(), newshapevec.end(), [&](      TableRow<2>& pt) { if (std::abs(pt.first) >= limitingHalfWidth_ ) pt.second.at(0) = extFactor_/bins; } );

    return newshapevec;

  }

  //============================================================================================================
  std::vector<double> ProtonPulseRandPDF::setSpectrum() {

    std::vector<double> spectrum;
    spectrum.reserve( times_.size() );

    if ( pulseEnum_ == ALLFLAT ) {
      spectrum = std::vector<double>( times_.size() , 1. );
    }
    else {

      for ( const auto& t : times_ )
        spectrum.push_back( pulseShape_.getValueAtKey(t).at(0)*acdipole_.getValueAtKey(t).at(0) );

    }
    return spectrum;
  }

  //============================================================================================================
  double ProtonPulseRandPDF::determineIntrinsicExt( const std::string& shapeTxtFile ) {
    Table<2> potShape = loadTable<2>( shapeTxtFile );
    auto const& shapevec = potShape.getShape( potShape(0,0), potShape( potShape.getNrows()-1, 0 ), tres_ );

    double num(0.), denom(0.);
    std::for_each( shapevec.begin(), shapevec.end(),
                   [&](const TableRow<2>& pt )
                   {
                     if ( std::abs(pt.first) >= limitingHalfWidth_ ) { num += pt.second.at(0);}
                     else denom += pt.second.at(0);
                   } );

    return num/denom;
  }

}
