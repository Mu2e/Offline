//
// Build a dictionary.
//
// $Id: classes.h,v 1.4 2011/06/11 01:50:47 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/11 01:50:47 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "art/Persistency/Common/Wrapper.h"

#include "Sandbox/inc/TransientProduct00Collection.hh"
#include "Sandbox/inc/TracerProduct.hh"
#include "Sandbox/inc/TracerProductCollection.hh"
#include "GeneralUtilities/inc/BarePointerCollection.hh"

class TrkRecoTrk{
 public:
  TrkRecoTrk(){}
  ~TrkRecoTrk(){}
};

template class art::Wrapper<mu2e::TransientProduct00Collection>;
template class art::Wrapper<mu2e::TracerProduct>;
template class art::Wrapper<std::vector<mu2e::TracerProduct> >;
template class art::Wrapper<mu2e::TracerProductCollection>;
template class art::Wrapper<mu2e::BarePointerCollection<TrkRecoTrk> >;
