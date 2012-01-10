void mu2e_gdml()
{
  TGeoManager *geom = TGeoManager::Import("mu2e.gdml");
  geom->Export("mu2e.C");
  TBrowser *b = new TBrowser();
  geom->GetTopVolume()->Draw();
}

