//
// A script to find the volumes that are on the STM beamline from z = 12000 mm to z = 45000 mm
//

#include "find_volume_at_point.C"

void find_stm_beamline_volumes(TString gdmlname = "mu2e_common.gdml") {

  TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); // default is kRootUnits which is cm, s, GeV; kG4Units are mm, ns, MeV
  TGeoManager *geom = TGeoManager::Import(gdmlname);

  double x = -3904;
  double y = 0;
  for (double z = 12000; z < 45000; z += 10) {
    find_volume_at_point(x, y, z, geom);
  }
}
