//
// Build a dictionary.
//
// $Id: classes.h,v 1.11 2012/04/27 05:37:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/27 05:37:32 $
//
// Original author Rob Kutschke
//

#include "DataProducts/inc/StrawIndex.hh"
#include "DataProducts/inc/FilterEfficiency.hh"

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include "cpp0x/array"
#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/EulerAngles.h"
#include <CLHEP/Geometry/Transform3D.h>

template class std::vector<CLHEP::Hep2Vector>;
template class std::vector<cet::map_vector_key>;
template class art::Wrapper<mu2e::FilterEfficiency>;

template class std::array<double,5>; // used in TubsParams
