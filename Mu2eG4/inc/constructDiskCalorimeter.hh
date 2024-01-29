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
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4PVPlacement.hh"



namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructDiskCalorimeter(VolumeInfo const& mother, SimpleConfig const& config);

  VolumeInfo caloBuildDisk      (const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildFrontPlate(const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildCase      (const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildBackPlate (const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildFEB       (const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildCrate     (const SimpleConfig& config, unsigned idisk);
  VolumeInfo caloBuildCable     (const SimpleConfig& config, unsigned idisk, VolumeInfo FEBvol);

  std::vector<G4TwoVector> caloExtrudedVertices(const std::vector<double>& stepsX, const std::vector<double>& stepsY, double delta=0.0);
  std::vector<G4double>    calcFEBPhiRange     (const DiskCalorimeter& cal);
  const G4LogicalVolume*   findCaloSolid       (const G4LogicalVolume* volume, G4String objectName, std::vector<const G4LogicalVolume*>& nodes);
  G4LogicalVolume*         caloLogical         (VolumeInfo volume, G4Material* mat, bool isVisible, const G4Color& color, bool isSolid, bool forceEdge);
  G4PVPlacement*           caloPlacement       (VolumeInfo& volume, const VolumeInfo& parent, G4RotationMatrix* rot, const G4ThreeVector& position,
                                                bool pMany, int copyNo, const SimpleConfig& config, bool doSurfaceCheck, int verbosityLevel);

}

#endif
