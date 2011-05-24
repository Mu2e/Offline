#ifndef RecoDataProducts_HoughCircleCollection_hh
#define RecoDataProducts_HoughCircleCollection_hh

//
// Define a type for a collection of Hough Circle Objects
//
// $Id: HoughCircleCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Peter Shanahan
//

#include <vector>

#include "RecoDataProducts/inc/HoughCircle.hh"

namespace mu2e {
   typedef std::vector<mu2e::HoughCircle> HoughCircleCollection;
}

#endif /* RecoDataProducts_HoughCircleCollection_hh */
