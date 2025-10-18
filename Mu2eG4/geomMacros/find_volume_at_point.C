//
// Useful functions to find a volume at a specific point
// see www.github.com/Mu2e/Offline/blob/main/Mu2eG4/geomMacros/README.md for important details
//

#include "clean_volume_names.C"

void find_volume_at_point(double x, double y, double z, TGeoManager* geom, bool output_csv = false, std::ostream& out = std::cout) { // can pass in your own geom
  TGeoNode * node = geom->FindNode(x, y, z);
  if (!node) out << "The position (" << x << ", " << y << ", " << z << ") is outside the world\n";
  else {
    if (output_csv) {
      out << x << "," << y << "," << z << "," << clean_volume_names(std::string(node->GetName())) << "," << node->GetVolume()->GetMaterial()->GetName() << "," << clean_path_names(std::string(geom->GetPath())) << "\n";
    }
    else {
      out << "The volume at position (" << x << ", " << y << ", " << z << ") is " << node->GetName() << "  (material = " << node->GetVolume()->GetMaterial()->GetName() << ", geometry path: " << geom->GetPath() << ")\n";
    }
  }
}

void find_volume_at_point(double x, double y, double z, TString gdmlname = "mu2e_common.gdml", bool output_csv = false, std::ostream& out = std::cout) {

  TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); // default is kRootUnits which is cm, s, GeV; kG4Units are mm, ns, MeV
  TGeoManager* geom = TGeoManager::Import(gdmlname);

  find_volume_at_point(x, y, z, geom, output_csv, out);
}

void find_volume_at_point_csv_header(std::ostream& out = std::cout) {
  out << "x,y,z,node,material,path\n";
}
