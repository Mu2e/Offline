// Explicit specializations for Table template class.  Includes
// function to interpolate tabulated data given a desired range of
// keys that are already given as a std::vector<double> object, or can
// be created based on a given key range and resolution.
//
// $Id: Table.cc,v 1.2 2014/04/25 17:44:13 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/04/25 17:44:13 $
//
// Original author: Kyle Knoepfel

#include "Mu2eUtilities/inc/Table.hh"
#include <numeric>

namespace mu2e {

  //=========================================================================
  template<> Value<1> Table<2>::getValueAtKey( const double key, const double binCorr ) const {

    // Use linear interpolation algorithm
    //
    // - note that in the case of pdfs, if the desired set of keys
    // does not match that of the stored table, it is the
    // responsibility of the user to ensure that the resulting pdf
    // integrates to 1.  This typically involves introducing a
    // bin-width correction factor to account for the difference in
    // key ranges and values.

    Value<1> val{{ 0. }};
    if      ( key == rawTable_.begin()->first ) val.at(0) = rawTable_.begin()->second.at(0);
    else if ( key >  rawTable_.begin()->first && key <= std::prev( rawTable_.end() )->first ) {
      auto const& it1 = this->getLowerBoundRow( key );
      auto const& it0 = std::prev( it1 );
      val.at(0)  = it0->second.at(0) + (it1->second.at(0) - it0->second.at(0))/(it1->first - it0->first)*(key-it0->first);
    }
    val.at(0) *= binCorr;
    return val;
  }

  //=========================================================================
  template<> void Table<2>::interpolateShape( const std::vector<double>& keys, const double binCorr ) {

    TableVec<2> shape;
    for ( const auto& k : keys )
      shape.emplace_back( k, getValueAtKey( k, binCorr ) );

    this->replaceShape( shape );

  }

  //=========================================================================
  template<> TableVec<2> Table<2>::getShape( const std::vector<double>& keys, const double binCorr ) const {

    TableVec<2> shape;
    for ( const auto& k : keys )
      shape.emplace_back( k, getValueAtKey( k, binCorr ) );

    return shape;
  }

  //=========================================================================
    template<> TableVec<2> Table<2>::getShape( const double keyLow, const double keyHigh, const double res, const double binCorr ) const {

    TableVec<2> shape;
    for ( double k(keyLow) ; k<=keyHigh ; k += res )
      shape.emplace_back( k, getValueAtKey( k, binCorr ) );

    return shape;
  }

  //=========================================================================
  template<> void Table<2>::renormalizeShape( const double norm ) {

    static auto add_to_integral = [](double integral, const TableRow<2>& it) {
      return integral + it.second.at(0);
    };

    const double integral = std::accumulate( rawTable_.begin(), rawTable_.end(), 0., add_to_integral );
    std::for_each( rawTable_.begin(), rawTable_.end(), [&](TableRow<2>& pt){ pt.second.at(0) *= norm/integral; } );

  }

} // end of namespace mu2e
