//
// Build a dictionary.
//
// $Id: classes.h,v 1.9 2012/04/05 18:14:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/05 18:14:35 $
//
// Original author Rob Kutschke
//

#include "DataProducts/inc/StrawIndex.hh"
#include "DataProducts/inc/FilterEfficiency.hh"

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/EulerAngles.h"

template class std::vector<cet::map_vector_key>;
template class art::Wrapper<mu2e::FilterEfficiency>;
