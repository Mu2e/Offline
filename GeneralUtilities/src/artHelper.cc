// Some helper functions for art
//
// $Id: artHelper.cc,v 1.2 2014/03/31 15:16:54 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/03/31 15:16:54 $
//
// Original author: Kyle Knoepfel

// Mu2e includes
#include "GeneralUtilities/inc/artHelper.hh"

namespace mu2e {

  std::vector<art::InputTag> artInputTagVector( const std::vector<std::string>& sv,
                                                const std::string str) {

    std::vector<art::InputTag> itv;
    itv.reserve( sv.size() );
    
    for ( auto const& ist : sv )
      itv.emplace_back( ist, str ); // implicit conversion to art::InputTag

    return itv;

  }

}
