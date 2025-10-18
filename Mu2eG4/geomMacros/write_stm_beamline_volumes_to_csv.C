//
// A script to find the volumes that are on the STM beamline from z = 12000 mm to z = 45000 mm
// and writes them out to a csv
//

#include "find_volume_at_point.C"

void write_stm_beamline_volumes_to_csv(TString gdmlname = "mu2e_common.gdml", TString csvname = "out.csv") {

  TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); // default is kRootUnits which is cm, s, GeV; kG4Units are mm, ns, MeV
  TGeoManager *geom = TGeoManager::Import(gdmlname);

  std::ofstream fout(csvname);
  bool output_csv = true;
  find_volume_at_point_csv_header(fout);

  double x = -3904;
  double y = 0;
  for (double z = 12000; z < 45000; z += 10) {
    find_volume_at_point(x, y, z, geom, output_csv, fout);
  }
}
