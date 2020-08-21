#ifndef Mu2eG4_constructTS_hh
#define Mu2eG4_constructTS_hh
//
// Free function to create  Transport Solenoid
//
//
// Original author KLG
//
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Material.hh"

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;
  class Beamline;

  void constructTS         ( VolumeInfo const& p, SimpleConfig const& c);

  void constructCryostat   ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCoils      ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCAs        ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructCollimators( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructDegrader   ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void constructPbarWindow ( VolumeInfo const& p, SimpleConfig const& c, Beamline const& bl);
  void addThermalShield    ( TransportSolenoid const& ts, VolumeInfo const& useAsParent, SimpleConfig const& c,
			     TransportSolenoid::TSRegion::enum_type TSRegion,
			     G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldAlMaterial,
			     double centerWallThickness);
  void addThermalShieldStraightSection( VolumeInfo const& useAsParent, std::vector<double> innerRadii,
					std::vector<double> outerRadii, double halfLength,
					G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldMidMaterial,
					CLHEP::Hep3Vector const& origin, std::string name);
  void addThermalShieldTorusSection( VolumeInfo const& useAsParent, std::vector<double> innerRadii,
				     std::vector<double> outerRadii, std::vector<double> torusParams,
				     G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldMidMaterial,
				     CLHEP::Hep3Vector const& origin, std::string name);

}

#endif /* Mu2eG4_constructTS_hh */
