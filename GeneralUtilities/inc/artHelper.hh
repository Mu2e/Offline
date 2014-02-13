#ifndef GeneralUtilities_artHelper_hh
#define GeneralUtilities_artHelper_hh

// Some helpers for interfacing art with Mu2e code
//
// $Id: artHelper.hh,v 1.1 2014/02/13 18:52:15 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/13 18:52:15 $
//
// Original author: Kyle Knoepfel

// C++ includes
#include <vector>

// Framework includes
#include "art/Utilities/InputTag.h"

namespace mu2e {

  // Convert vector of std::string's to vector of art::InputTag's
  std::vector<art::InputTag> artInputTagVector( const std::vector<std::string>& sv ); 

}

#endif
