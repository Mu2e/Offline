#ifndef GeneralUtilities_artHelper_hh
#define GeneralUtilities_artHelper_hh

// Some helpers for interfacing art with Mu2e code
//
// $Id: artHelper.hh,v 1.2 2014/03/31 15:16:54 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/03/31 15:16:54 $
//
// Original author: Kyle Knoepfel

// C++ includes
#include <vector>

// Framework includes
#include "art/Utilities/InputTag.h"

namespace mu2e {

  // Convert vector of std::string's to vector of art::InputTag's
  std::vector<art::InputTag> artInputTagVector( const std::vector<std::string>& sv,
                                                const std::string str = "" ); 

}

#endif
