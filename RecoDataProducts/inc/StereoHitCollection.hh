#ifndef RecoDataProducts_StereoHitCollection_hh
#define RecoDataProducts_StereoHitCollection_hh

//
// Define a type for a collection of StereoHit objects.
//
// $Id: StereoHitCollection.hh,v 1.1 2013/01/26 18:18:44 brownd Exp $
// $Author: brownd $
// $Date: 2013/01/26 18:18:44 $
//
// Original author David Brown
//

#include <vector>

#include "RecoDataProducts/inc/StereoHit.hh"

namespace mu2e {
   typedef std::vector<mu2e::StereoHit> StereoHitCollection;
}

#endif /* RecoDataProducts_StereoHitCollection_hh */
