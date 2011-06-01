//
// Build a dictionary.
//
// $Id: classes.h,v 1.4 2011/06/01 14:57:48 greenc Exp $
// $Author: greenc $
// $Date: 2011/06/01 14:57:48 $
//
// Original author Rob Kutschke
//

#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "TrackerGeom/inc/StrawId.hh"
#include "TrackerGeom/inc/StrawIndex.hh"
#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Provenance/ProductID.h"
#include "cetlib/map_vector.h"
#include <vector>

template class art::Wrapper<mu2e::DPIndexVectorCollection>;

template class std::vector<cet::map_vector_key>;
template class std::vector<mu2e::DPIndex>;
template class std::vector<mu2e::DPIndexVector>;
