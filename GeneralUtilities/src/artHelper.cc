// Some helper functions for art
//
// $Id: artHelper.cc,v 1.1 2014/02/13 18:52:15 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/13 18:52:15 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "GeneralUtilities/inc/artHelper.hh"

namespace mu2e {

  std::vector<art::InputTag> artInputTagVector( const std::vector<std::string>& sv ) {

    std::vector<art::InputTag> itv;
    itv.reserve( sv.size() );
    
    for ( auto const& ist : sv )
      itv.push_back( ist ); // implicit conversion to art::InputTag

    return itv;

  }

}
