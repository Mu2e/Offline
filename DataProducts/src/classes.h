//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#include "canvas/Persistency/Common/Wrapper.h"
#include "cetlib/map_vector.h"
#include "boost/array.hpp"
#include <vector>
#include <utility>
#include <TString.h>

// PDG
#include "Offline/DataProducts/inc/PDGCode.hh"

// straws
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/DataProducts/inc/PanelId.hh"
#include "Offline/DataProducts/inc/LayerId.hh"
#include "Offline/DataProducts/inc/PlaneId.hh"

// tracker
#include "Offline/DataProducts/inc/Helicity.hh"

// calorimeter
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"

// CRV
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

// ExtMon
#include "Offline/DataProducts/inc/ExtMonFNALModuleId.hh"
#include "Offline/DataProducts/inc/ExtMonFNALChipId.hh"
#include "Offline/DataProducts/inc/ExtMonFNALPixelId.hh"

// General
#include "Offline/DataProducts/inc/IndexMap.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

// trigger
#include "Offline/DataProducts/inc/EventWindowMarker.hh"

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

// STM
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/DataProducts/inc/STMTestBeamEventInfo.hh"

// General
#include "Offline/DataProducts/inc/SurfaceId.hh"
