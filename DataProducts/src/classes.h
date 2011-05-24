//
// Build a dictionary.
//
// $Id: classes.h,v 1.3 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "TrackerGeom/inc/StrawId.hh"
#include "TrackerGeom/inc/StrawIndex.hh"
#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Provenance/ProductID.h"
#include "cetlib/map_vector.h"
#include <string>
#include <vector>

template class art::Wrapper<mu2e::DPIndexVectorCollection>;

template class std::vector<cet::map_vector_key>;
template class std::vector<mu2e::DPIndex>;
template class std::vector<mu2e::DPIndexVector>;
