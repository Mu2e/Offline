#ifndef ToyDP_HoughCircleCollection_hh
#define ToyDP_HoughCircleCollection_hh

//
// Define a type for a collection of ToyHits.
//
// $Id: HoughCircleCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Peter Shanahan
//

#include <vector>

#include "ToyDP/inc/HoughCircle.hh"

namespace mu2e {
   typedef std::vector<mu2e::HoughCircle> HoughCircleCollection;
}

#endif /* ToyDP_HoughCircleCollection_hh */
