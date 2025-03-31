//
// Root macro to draw an existing volume with its 3D projections and vertices
//
// Inputs:
// -------
//
//   - The volume name (TString)
//
//   - The gdml file prefix
//
// Key features:
// -------------
//
//    - Box origin:
//        Get the origin of the wanted volume in local coordinates
//        Loop on the nodes of geometry top volume to get the local to global
//        coordinates transformation matrix for the wanted volume
//        Transform the local coordinates to global coordinates (in mm)
//
//    - Volumes Visualization:
//        Make volume visible in the ROOT OpenGL viewer.
//
// Output:
// -------
//
//    - Opens a TBrowser for interactive exploration.
//
// How to run:
// -----------
//
//    root -l draw_volume.C\(\"[volume_name]\",\"[gdml_file_prefix]\"\)
//
//    Example:
//    root -l draw_volume.C\(\"SRetainingWallFoot\",\"mu2e_40\"\)
//

#include <iostream>
#include <algorithm>
#include "TGeoVolume.h"
#include "TBrowser.h"

int draw_volume(
    const TString& volName = "SRetainingWallFoot",
    const TString& gdmlFilePrefix = "mu2e_40") {

  // ********** IMPORT GEOMETRY **********
  TString gdmlFile = gdmlFilePrefix + ".gdml";
  if (gSystem->AccessPathName(gdmlFile.Data())) {
    std::cerr << "Error: GDML file not found: " << gdmlFile << std::endl;
    return -1;
  }

  TGeoManager* geom = TGeoManager::Import(gdmlFile);
  if (!geom) {
    std::cerr << "Error: Unable to load geometry from GDML file: " << gdmlFile << std::endl;
    return -2;
  }

  // ********** CHECK IF THE VOLUME EXISTS **********
  TGeoVolume* volume = geom->FindVolumeFast(volName.Data());
  if (!volume) {
    std::cerr << "Error: Volume not found: " << volName << std::endl;
    return 1;
  }

  // ********** GET VOLUME PARAMETERS **********
  TGeoBBox* bbox = dynamic_cast<TGeoBBox*>(volume->GetShape());
  if (!bbox) {
    std::cerr << "Error: Unable to retrieve bounding box for volume: " << volName << std::endl;
    return 2;
  }

  // Local bounding box information
  const double* localOrigin = bbox->GetOrigin();  // Local origin in cm
  double localHalfLengths[3] = {bbox->GetDX(), bbox->GetDY(), bbox->GetDZ()};
  double localEdges[3] = {
      localOrigin[0] + localHalfLengths[0],
      localOrigin[1] + localHalfLengths[1],
      localOrigin[2] + localHalfLengths[2]};

  // ********** FIND NODE IN GEOMETRY TREE **********
  TGeoIterator iterator(geom->GetTopVolume());
  TGeoNode* currentNode = nullptr;
  bool nodeFound = false;

  while ((currentNode = iterator())) {
    if (volName == currentNode->GetVolume()->GetName()) {
      nodeFound = true;
      break;
    }
  }

  if (!nodeFound) {
    std::cerr << "Error: Node named " << volName << " not found in the geometry tree." << std::endl;
    return -3;
  }

  // ********** TRANSFORM COORDINATES TO GLOBAL **********
  const TGeoMatrix* localToGlobal = iterator.GetCurrentMatrix();
  if (!localToGlobal) {
    std::cerr << "Error: Unable to retrieve local-to-global transformation matrix." << std::endl;
    return -4;
  }

  double globalOrigin[3];
  double globalEdges[3];
  double globalHalfLengths[3];

  // Transform local to global coordinates
  localToGlobal->LocalToMaster(localOrigin, globalOrigin);
  localToGlobal->LocalToMaster(localEdges, globalEdges);

  // Calculate global half lengths
  for (int i = 0; i < 3; ++i) {
    globalHalfLengths[i] = globalEdges[i] - globalOrigin[i];
  }

  // ********** PRINT PARAMETERS **********
  std::cout << "Volume Name: " << currentNode->GetVolume()->GetName() << std::endl;
  std::cout << "Global Position (mm): (" << globalOrigin[0] * 10 << ", "
            << globalOrigin[1] * 10 << ", " << globalOrigin[2] * 10 << ")" << std::endl;
  std::cout << "Local Bounding Box Half Lengths (mm): ("
            << localHalfLengths[0] * 10 << ", " << localHalfLengths[1] * 10
            << ", " << localHalfLengths[2] * 10 << ")" << std::endl;
  std::cout << "Global Bounding Box Half Lengths (mm): ("
            << globalHalfLengths[0] * 10 << ", " << globalHalfLengths[1] * 10
            << ", " << globalHalfLengths[2] * 10 << ")" << std::endl;

  // ********** DRAW VOLUME **********
  // Hide all volumes
  std::for_each(geom->GetListOfVolumes()->begin(),
                geom->GetListOfVolumes()->end(),
                [](TObject* obj) {
                  static_cast<TGeoVolume*>(obj)->SetVisibility(0);
                });

  // Highlight the selected volume
  volume->SetVisibility(1);
  volume->SetLineColor(kRed);

  // Open a browser and draw the geometry
  TBrowser* browser = new TBrowser();
  geom->GetVolume("HallAir")->Draw("ogl");

  return 0;
}
