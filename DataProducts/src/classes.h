//
// Build a dictionary.
//
// $Id: classes.h,v 1.18 2013/07/01 06:52:13 gandr Exp $
// $Author: gandr $
// $Date: 2013/07/01 06:52:13 $
//
// Original author Rob Kutschke
//

#include "DataProducts/inc/StrawIndex.hh"
#include "DataProducts/inc/FilterEfficiency.hh"
#include "DataProducts/inc/ExtMonFNALSensorId.hh"
#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "DataProducts/inc/ExtMonFNALChipId.hh"
#include "DataProducts/inc/ExtMonFNALPixelId.hh"

#include "art/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include "boost/array.hpp"
#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <CLHEP/Geometry/Transform3D.h>

template class std::vector<CLHEP::Hep2Vector>;
template class std::vector<cet::map_vector_key>;
template class art::Wrapper<mu2e::FilterEfficiency>;

template class boost::array<double,5>; // used in TubsParams
