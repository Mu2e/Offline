//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2012/03/05 20:14:15 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/05 20:14:15 $
//
// Original author Rob Kutschke
//

#include "DataProducts/inc/StrawIndex.hh"
#include "DataProducts/inc/FilterEfficiency.hh"

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include <vector>

template class std::vector<cet::map_vector_key>;
template class art::Wrapper<mu2e::FilterEfficiency>;
