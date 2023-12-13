#ifndef Mu2eG4_constructDiskCalorimeter_hh
#define Mu2eG4_constructDiskCalorimeter_hh
//
// Free function to create the disk calorimeter.
//
//
// Original author Rob Kutschke
//
// Notes:
// 1) Arguments are:
//    1 - pointer to the mother logical volume.
//    2 - geometry file

// Mu2e includes.
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4TwoVector.hh"



namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructDiskCalorimeter(VolumeInfo const& mother, SimpleConfig const& config);

  G4LogicalVolume* caloBuildLogical(G4VSolid* solid, G4Material* mat, const G4String& name, bool isVisible, const G4Color&  color, bool isSolid, bool forceEdge);
  G4LogicalVolume* caloBuildFrontPlate(const SimpleConfig& config,MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildDisk(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildBackPlate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildCrate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal);
  G4LogicalVolume* caloBuildFEB(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal);
  std::vector<G4TwoVector> caloExtrudedVertices(const std::vector<double>& stepsX, const std::vector<double>& stepsY, double delta=0.0);


}

#endif /* Mu2eG4_constructDiskCalorimeter_hh */
