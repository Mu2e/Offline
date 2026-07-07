//
// A macro to find the volume at a given point
// see www.github.com/Mu2e/Offline/blob/main/Mu2eG4/geomMacros/README.md for important details
//
void where_is_volume(TString search_name, TString gdmlname = "mu2e_common.gdml", bool verbose = false) { // don't need to know the full real name, will look for volumes that contain this string
  TGeoManager::SetDefaultUnits(TGeoManager::kG4Units); // default is kRootUnits which is cm, s, GeV; kG4Units are mm, ns, MeV

  TGeoManager *geom = TGeoManager::Import(gdmlname);

  TGeoIterator next(geom->GetTopVolume());
  TGeoNode*current; TGeoVolume *vol;
  TString path;

  while ((current=next())) {
    vol = current->GetVolume();
    next.GetPath(path);
    auto global = next.GetCurrentMatrix();
    // if you want to see where is the center (origin) of this object in TOP coordinates:
    const double *local_bbox_orig = ((TGeoBBox*)vol->GetShape())->GetOrigin();
    double master_orig[3];
    global->LocalToMaster(local_bbox_orig, master_orig);
    TString vol_name = vol->GetName();
    if (vol_name.Contains(search_name)) {
      std::cout << vol_name << " is at global position: (" << master_orig[0] <<   ", " << master_orig[1] << ", " << master_orig[2] << ")\n";
      if (verbose) {
        std::cout << "\t geometry path: " << path << "\n";
        std::cout << "\t material: " << vol->GetMaterial()->GetName() << "\n";
        std::cout << "\t shape: ";vol->GetShape()->Dump();
        std::cout << std::endl;
      }
    }
  }
}
