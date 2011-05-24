//
// Build a dictionary.
//
// $Id: classes.h,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Provenance/ProductID.h"

#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "GeneralUtilities/inc/MapVectorKey.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "TrackerGeom/inc/StrawId.hh"
#include "TrackerGeom/inc/StrawIndex.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/LorentzVector.h"

template class art::Wrapper<mu2e::DPIndexVectorCollection>;

template class std::vector<MapVectorKey>;
template class std::vector<mu2e::DPIndex>;
template class std::vector<mu2e::DPIndexVector>;
template class std::vector<CLHEP::Hep3Vector>;

