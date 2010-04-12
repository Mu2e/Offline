#ifndef ToyDP_HoughCircleCollection_hh
#define ToyDP_HoughCircleCollection_hh

//
// Define a type for a collection of ToyHits.
//
// $Id: HoughCircleCollection.hh,v 1.1 2010/04/12 18:12:28 shanahan Exp $
// $Author: shanahan $
// $Date: 2010/04/12 18:12:28 $
//
// Original author Peter Shanahan
//

#include <vector>

#include "ToyDP/inc/HoughCircle.hh"

namespace mu2e {
   typedef std::vector<mu2e::HoughCircle> HoughCircleCollection;
}

#endif //ToyDP_HoughCircleCollection_hh
