void overlapCheck(TString const fname="mu2e.gdml",
                  Double_t res=0.001)
{
  TGeoManager::Import(fname);
  gGeoManager->CheckOverlaps(res);
  gGeoManager->PrintOverlaps();
  return;
}
