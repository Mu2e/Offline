//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2011/06/05 17:29:11 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/05 17:29:11 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "Sandbox/inc/TransientProduct00Collection.hh"
#include "Sandbox/inc/TracerProduct.hh"

template class art::Wrapper<mu2e::TransientProduct00Collection>;
template class art::Wrapper<mu2e::TracerProduct>;
template class art::Wrapper<std::vector<mu2e::TracerProduct> >;
