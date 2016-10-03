// A carbon plane with several modules attatched.
//
// Evan Schiewe, 2013

#ifndef EXTMONFNALPLANE_HH
#define EXTMONFNALPLANE_HH

#include "canvas/Persistency/Common/Wrapper.h"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "DataProducts/inc/ExtMonFNALModuleDenseId.hh"
#include "DataProducts/inc/ExtMonFNALPixelId.hh"
#include "DataProducts/inc/ExtMonFNALPlaneId.hh"

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }
  namespace ExtMonFNAL { class ExtMonMaker; }

  class ExtMonFNALPlane {
  public:
    // the plane itself
    std::vector<double> halfSize() const { return halfSize_; }

    // contains several modules
    const std::vector<double>& module_zoffset() const { return m_module_zoffset; }
    const std::vector<double>& module_xoffset() const { return m_module_xoffset; }
    const std::vector<double>& module_yoffset() const { return m_module_yoffset; }
    const std::vector<double>& module_rotation() const { return m_module_rotation; }

    unsigned nModules() const { return m_module_zoffset.size(); }

    //Not implemented - comment out to satisfy ROOT dictionary maker.
    //ExtMonFNALPixelId findPixel(ExtMonFNALPlaneId pid, double xPlane, double yPlane) const;

    const ExtMonFNALModule& module() const { return module_; }

    CLHEP::Hep3Vector planeCoordinates(const ExtMonFNALPixelId& pix) const;
    const std::string& planeMaterial() const { return m_planeMaterial; }

    ExtMonFNALPlane(const ExtMonFNALModule& module, const std::vector<double>& hs)
      : module_(module), halfSize_(hs)

    {}
    // Required by genreflex persistency
    ExtMonFNALPlane() {}
    //---------------------------------------------------------------------
  private:
    ExtMonFNALModule module_;
    std::vector<double> halfSize_;

    // Plane material
    std::string m_planeMaterial;

    // Module center positions
    std::vector<double> m_module_zoffset;
    std::vector<double> m_module_xoffset;
    std::vector<double> m_module_yoffset;
    std::vector<double> m_module_rotation;

    template<class T> friend class art::Wrapper;

    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNAL::ExtMonMaker;
  };


} // namespace mu2e

#endif /*EXTMONFNALPLANE_HH*/
