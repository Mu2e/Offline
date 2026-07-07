// A macro to draw detector from GDML in ROOT
// Author : Zhengyun You
//
// -------------------------------------------------------------------------------------------------------------
//
// 1. In left column, select Default->Master Volume-> In tree of volumes, right click the one to draw -> Draw
// 2. Tick on/off to show/hide the volumes
// 3. For 3D OpenGL view, in "View" tab -> "View with" -> "OpenGL"
// 4. To open a box on detector to see inside in GL viewer, select "Clipping" -> "Box",
//    select "Edit in Viewer" or change values and "Apply".
// 5. The color is not consistent with that defined in G4, to change colors, edit the output "mu2e.C" and run it.
//
// -------------------------------------------------------------------------------------------------------------

TBrowser *b = 0;
void mu2e_gdml()
{
  TGeoManager *geom = TGeoManager::Import("mu2e_common.gdml");
  geom->Export("mu2e.C");
  b = new TBrowser();
  geom->GetTopVolume()->Draw();
}
