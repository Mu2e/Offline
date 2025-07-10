void overlapCheck(TString const fname="mu2e.gdml",
                  Double_t res=1.e-12)
{
  TGeoManager::Import(fname);
  gGeoManager->CheckOverlaps(res);
  gGeoManager->PrintOverlaps();
  return;
}
