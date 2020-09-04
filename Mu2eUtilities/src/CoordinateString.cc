//
// Original author: Kyle Knoepfel

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

#include "GeneralUtilities/inc/splitLine.hh"

#include "Mu2eUtilities/inc/CoordinateString.hh"

namespace mu2e {

  //=========================================================================
  CoordinateString::CoordinateString( const std::string& coordStr )
    : coordStr_( coordStr             )
    , coordStd_( readCoordinatesStd( coordStr_ ) )
    , coord_   ( calcCoordinates( coordStd_ )    )
  {}


  //=========================================================================
  std::array<CoordinateString::FtInchPair,2> CoordinateString::readCoordinatesStd( const std::string& coordStr ) {

    std::vector<std::string> parts;
    splitLine( coordStr, ",", parts );

    FtInchPair xStd = makeFtInchPair( parts.at(0) );
    FtInchPair yStd = makeFtInchPair( parts.at(1) );

    return std::array<FtInchPair,2>{{xStd, yStd}};
  }



  //=========================================================================
  std::array<double,2> CoordinateString::calcCoordinates( const std::array<FtInchPair,2>& coordStd ) {

    const double x = convert2mm( coordStd.at(0) );
    const double y = convert2mm( coordStd.at(1) );

    return std::array<double,2>{{x,y}};
  }

  //=========================================================================
  CoordinateString::FtInchPair CoordinateString::makeFtInchPair( const std::string& stringToParse ) {
    std::vector<std::string> stringPairing;
    splitLine( stringToParse, ":", stringPairing );

    const int    ft = ( stringPairing.size() > 0 ) ? std::atoi( stringPairing.at(0).c_str() ) : 0;
    const double in = ( stringPairing.size() > 1 ) ? std::atof( stringPairing.at(1).c_str() ) : 0.;

    return FtInchPair( ft, in );
  }

  //=========================================================================
  double CoordinateString::convert2mm( const FtInchPair& ftInchPair ) {
    double inches = ftInchPair.first*12;
    inches += ftInchPair.second;
    return inches*25.4;
  }

  //=========================================================================
  void CoordinateString::print() const {
    std::cout << " Location: ( " << coord_.at(0) << "," << coord_.at(1) << " ) or [ "
              << coordStd_.at(0).first << " ft " << coordStd_.at(0).second << " in , "
              << coordStd_.at(1).first << " ft " << coordStd_.at(1).second << " in ] " << std::endl;
  }

} // end of namespace mu2e


