//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#include "canvas/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include "boost/array.hpp"
#include <vector>

// straws
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include "DataProducts/inc/PanelId.hh"
#include "DataProducts/inc/LayerId.hh"
#include "DataProducts/inc/PlaneId.hh"

// tracker
#include "DataProducts/inc/Helicity.hh"

// CRV
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

// ExtMon
#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "DataProducts/inc/ExtMonFNALChipId.hh"
#include "DataProducts/inc/ExtMonFNALPixelId.hh"

// General
#include "DataProducts/inc/IndexMap.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

// trigger
#include "DataProducts/inc/EventWindowMarker.hh"

// CLHEP
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/EulerAngles.h"
#include <CLHEP/Geometry/Transform3D.h>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
