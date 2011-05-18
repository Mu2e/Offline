#ifndef ToyDP_HoughCircleCollection_hh
#define ToyDP_HoughCircleCollection_hh

//
// Define a type for a collection of Hough Circle Objects
//
// $Id: HoughCircleCollection.hh,v 1.3 2011/05/18 18:12:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 18:12:40 $
//
// Original author Peter Shanahan
//

#include <vector>

#include "ToyDP/inc/HoughCircle.hh"

namespace mu2e {
   typedef std::vector<mu2e::HoughCircle> HoughCircleCollection;
}

#endif /* ToyDP_HoughCircleCollection_hh */
