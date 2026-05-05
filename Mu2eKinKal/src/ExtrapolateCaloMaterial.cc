#include "Offline/Mu2eKinKal/inc/ExtrapolateCaloMaterial.hh"
#include "Offline/KinKalGeom/inc/Calo.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <iostream>
#include <limits>

namespace mu2e {
  using KinKal::VEC3;
  using KKGeom::Calo;

  ExtrapolateCaloMaterial::ExtrapolateCaloMaterial(double maxdt, double dptol, double intertol,
                                                   Calo const& calo, int depth, int debug) :
    maxDt_(maxdt), dptol_(dptol), intertol_(intertol), fp_z_(std::numeric_limits<double>::max()),
    inter_(), surftype_(Unknown), intersection_found_(false), fpann_(nullptr), debug_(debug)
  {
    // Get disk 0 or 1 based on depth
    auto const& diskref = (depth == 0) ? calo.EMC_Disk_0_Front() : calo.EMC_Disk_1_Front();

    // Z position of front panel
    // Front panel dimensions from geometry:
    // FPInnerRadius = 336, FPOuterRadius = 680
    // Total FP thickness: foam (21.75) + carbon (3.0) = 24.75 mm
    double fp_inner_r = 336.0;
    double fp_outer_r = 680.0;
    double disk_z = diskref.center().Z();

    // Front panel surfaces: positioned ~50mm before disk (FOA position)
    // This should be read from geometry service TODO
    fp_z_ = disk_z - 50.0; // approximate position

    // Create annulus surface for front panel (combined foam + carbon layers)
    // Annulus constructor: Annulus(norm, udir, center, innerrad, outerrad)
    VEC3 norm(0, 0, 1);      // normal to disk (Z-direction)
    VEC3 udir(1, 0, 0);      // u-direction along disk (radial)
    VEC3 center_fp(0, 0, fp_z_);
    fpann_ = std::make_shared<KinKal::Annulus>(norm, udir, center_fp,
                                              fp_inner_r, fp_outer_r);

    if(debug_ == -300) {
      std::cout << "\n=== ExtrapolateCaloMaterial Initialization ==="<< std::endl;
      std::cout << "  Disk index: " << depth << std::endl;
      std::cout << "  Disk center Z: " << disk_z << " mm" << std::endl;
      std::cout << "  Front panel Z: " << fp_z_ << " mm (offset: " << (disk_z - fp_z_) << " mm upstream)" << std::endl;
      std::cout << "  Annulus geometry:" << std::endl;
      std::cout << "    Inner radius: " << fp_inner_r << " mm" << std::endl;
      std::cout << "    Outer radius: " << fp_outer_r << " mm" << std::endl;
      std::cout << "    Normal: (" << norm.X() << ", " << norm.Y() << ", " << norm.Z() << ")" << std::endl;
      std::cout << "    U-direction: (" << udir.X() << ", " << udir.Y() << ", " << udir.Z() << ")" << std::endl;
      std::cout << "  Material thickness: 21.75 mm foam + 3.0 mm carbon = 24.75 mm total" << std::endl;
      std::cout << "  Max extrapolation time: " << maxDt_ << " (unbounded if < 0)" << std::endl;
      std::cout << "  Intersection tolerance: " << intertol_ << " mm" << std::endl;
      std::cout << "=== Initialization Complete ===\n" << std::endl;
    }
  }

}
