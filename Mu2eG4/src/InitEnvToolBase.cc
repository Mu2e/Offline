///////////////////////////////////////////////////////////////////////////////
// base class for constructing non-Mu2e geometries to be studied using 
// Mu2e simulation framework
///////////////////////////////////////////////////////////////////////////////

#include "Mu2eG4/inc/InitEnvToolBase.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
//  geometry construction hook
//-----------------------------------------------------------------------------
  int InitEnvToolBase::construct(VolumeInfo const & ParentVInfo, SimpleConfig const& Config) {
    _name = "unknown";
    return 0;
  }

}
