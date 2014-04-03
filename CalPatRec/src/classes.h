//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2014/04/03 23:18:57 murat Exp $
// $Author: murat $
// $Date: 2014/04/03 23:18:57 $
//
// Original author Rob Kutschke
//

#include <vector>
//#include <array>

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"

#include "CalPatRec/inc/AlgorithmID.hh"
#include "CalPatRec/inc/AlgorithmIDCollection.hh"

// Cannot use the typedefs in here - not sure why.

template class art::Ptr<mu2e::AlgorithmID>;
template class std::vector<art::Ptr<mu2e::AlgorithmID> >;
template class art::Wrapper<mu2e::AlgorithmIDCollection>;
