//
// Build a dictionary.
//
// $Id: classes.h,v 1.5 2011/06/11 02:27:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 02:27:32 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "Sandbox/inc/TransientProduct00Collection.hh"
#include "Sandbox/inc/TracerProduct.hh"
#include "Sandbox/inc/TracerProductCollection.hh"

template class art::Wrapper<mu2e::TransientProduct00Collection>;
template class art::Wrapper<mu2e::TracerProduct>;
template class art::Wrapper<std::vector<mu2e::TracerProduct> >;
template class art::Wrapper<mu2e::TracerProductCollection>;

