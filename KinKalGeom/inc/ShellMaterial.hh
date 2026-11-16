//
// Class modeling a thin material as a surface, material, and thickness. Can be used to create a ShellXing in KKTrack
// Original author: David Brown (LBNL) 5/26
//
#ifndef KinKalGeom_ShellMaterial_hh
#define KinKalGeom_ShellMaterial_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/MatEnv/DetMaterial.hh"

namespace mu2e {
  class ShellMaterial {
    using SurfacePtr = std::shared_ptr<KinKal::Surface>;
    public:
    private:

    SurfacePtr surf_; // surface for this materia


  };
}
