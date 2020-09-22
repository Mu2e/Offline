#ifndef Mu2eUtilities_CoordinateString_hh
#define Mu2eUtilities_CoordinateString_hh
//
//
// Original author: Kyle Knoepfel

//
// Class for converting an (x,y) pair from standard units (formatted
// as a string), to metric units (mm).  The construction of the object looks like:
//
//     CoordinateString("xft:xin,yft:yin"), 
//
// where xft, yft are integers expressed within the string, and xin,
// yin are doubles expressed within the string.  Acceptable examples
// include:
//
//     CoordinateString("-4:3,5:6") // parsed as x = -4' 3" | y =  5'  6"
//     CoordinateString(":3,4")     // parsed as x =  0' 3" | y =  4'  0"
//
// Note the following example:
//
//     CoordinateString("2:,-:5")   // parsed as x =  2' 0" | y = -0' +5" 
//
// you probably don't want this, instead you'll have to do something like:
//
//     CoordinateString("2:,:-5")   // parsed as x =  2' 0" | y =  0' -5"
//
// Also note that the number preceding the : is converted to an
// integer, whereas the number afterward is cast into a double so that
// floating-point inch values are allowed.  Note that a double value
// for feet does not make sense.  For example
//
//     CoordinateString("2.5:2,-3:2") // parsed as x = 2' 2" | y = -3' 2"
//
// No free function currently exists for adding coordinates, but it
// would be straightforward to do if needed.

// C++ includes
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>

// Framework includes
#include "cetlib_except/exception.h"

namespace mu2e {
  
  class CoordinateString {
    
  public:
    
    typedef std::pair <int,double> FtInchPair;
    
    // Constructors
    explicit CoordinateString( const std::string& coordStr ); 
    
    // Accessors
    const std::array<FtInchPair,2>& getPairStd( ) const { return coordStd_; }
    const std::array<double,2> & getPair() const { return coord_; }
    
    // Metric (mm)
    double get(std::size_t i) const { return coord_.at(i); }
    
    double x() const { return get(0); }
    double y() const { return get(1); }
    
    void print() const;

    static std::array<FtInchPair,2> readCoordinatesStd( const std::string& coordStr );
    static std::array<double,2>     calcCoordinates   ( const std::array<FtInchPair,2>& coordStd );

    static FtInchPair makeFtInchPair( const std::string& stringToParse );
    static double     convert2mm    ( const FtInchPair& ftInchPair );

  private:
    
    std::string coordStr_;
    std::array<FtInchPair,2> coordStd_;
    std::array<double,2> coord_;
    
  };
 
} // end of namespace mu2e

#endif /* Mu2eUtilities_CoordinateString_hh */
