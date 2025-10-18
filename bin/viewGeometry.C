void viewGeometry(const std::string& gdmlfile="mu2e_common.gdml")
{
  TGeoManager *geom = TGeoManager::Import(gdmlfile.c_str());
  geom->GetVolume("HallAir")->Draw("ogl");
}
