#ifndef TrackerGeom_PanelEB_hh
#define TrackerGeom_PanelEB_hh
//
// Holds information about the tracker electronics board (used in G4)
//
#include "GeomPrimitives/inc/TubsParams.hh"
#include <string>

namespace mu2e {
  struct PanelEB{
    PanelEB() {}
    PanelEB(TubsParams const& key, TubsParams const& shield, std::string const& keymat, std::string const& shieldmat, double rotation) :
      _EBKey(key), _EBKeyShield(shield), _EBKeyMaterial(keymat), _EBKeyShieldMaterial(shieldmat), _EBKeyPhiExtraRotation(rotation) {}
    TubsParams  const& getEBKeyParams()           const { return _EBKey; }
    TubsParams  const& getEBKeyShieldParams()     const { return _EBKeyShield; }
    std::string const& getEBKeyMaterial()         const { return _EBKeyMaterial; }
    std::string const& getEBKeyShieldMaterial()   const { return _EBKeyShieldMaterial; }
    double             getEBKeyPhiExtraRotation() const { return _EBKeyPhiExtraRotation; }
    //
    TubsParams  _EBKey;
    TubsParams  _EBKeyShield;
    std::string _EBKeyMaterial;
    std::string _EBKeyShieldMaterial;
    double      _EBKeyPhiExtraRotation;
  };
}
#endif
