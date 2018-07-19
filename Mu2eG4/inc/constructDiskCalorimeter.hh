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
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "G4String.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4Color.hh"



namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructDiskCalorimeter(VolumeInfo const& mother, SimpleConfig const& config);

  G4LogicalVolume* caloBuildLogical(G4VSolid* solid, G4Material* mat, const G4String& name, bool isVisible, const G4Color&  color, bool isSolid, bool forceEdge);
  G4LogicalVolume* caloBuildFrontPlate(const SimpleConfig& config,MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildDisk(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildBackPlate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal, int idisk);
  G4LogicalVolume* caloBuildCrate(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal);
  G4LogicalVolume* caloBuildFEB(const SimpleConfig& config, MaterialFinder& materialFinder, const DiskCalorimeter& cal);


}

#endif /* Mu2eG4_constructDiskCalorimeter_hh */
