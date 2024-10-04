void find_volume_at_point(double x, double y, double z, TGeoManager* geom) { // can pass in your own geom
  TGeoNode * node = geom->FindNode(x, y, z);
  if (!node) std::cout << "The position (" << x << ", " << y << ", " << z << ") is outside the world\n";
  else {
    std::cout << "The volume at position (" << x << ", " << y << ", " << z << ") is " << node->GetName() << "  (material = " << node->GetVolume()->GetMaterial()->GetName() << ", geometry path: " << geom->GetPath() << ")\n";
  }
}

void find_volume_at_point(double x, double y, double z, TString gdmlname = "mu2e.gdml") {

  TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); // default is kRootUnits which is cm, s, GeV; kG4Units are mm, ns, MeV
  TGeoManager* geom = TGeoManager::Import(gdmlname);

  find_volume_at_point(x, y, z, geom);
}
