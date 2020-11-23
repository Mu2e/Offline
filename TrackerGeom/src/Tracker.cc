//
// Geometry and identifier info about an Tracker.
//
//
//
// Original author Rob Kutschke
//

#include <utility>
#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {

  Tracker::Tracker(StrawCollection const& straws, StrawProperties const& sprops, const char* name) :
  _name(name), _strawprops(sprops), _strawindex{}, _straws(straws) {
    // create the fast lookup map
    for(uint16_t plane=0; plane < StrawId::_nplanes; plane++){
      for(uint16_t panel = 0;panel < StrawId::_npanels; panel++){
	for(uint16_t straw = 0; straw < StrawId::_nstraws; straw++){
	  StrawId sid(plane,panel,straw);
	  _strawindex[sid.asUint16()] = sid.uniqueStraw();
	}
      }
    }
    // build the panels from the straws
    for(uint16_t plane=0; plane < StrawId::_nplanes; plane++){
      for(uint16_t panel = 0;panel < StrawId::_npanels; panel++){
	StrawId sid(plane,panel,0);
	_panels[sid.uniquePanel()] = Panel(sid, _straws);
      }
    }
    // build the planes from the panels
    for(uint16_t plane=0; plane < StrawId::_nplanes; plane++){
      StrawId sid(plane,0,0);
      _planes[sid.plane()] = Plane(sid, _panels);
    }

  }

// the following copies the core tracker plus all the elements built outside this class by TrackerMaker: this design needs to be refactored FIXME
  Tracker::Tracker(const Tracker& other) : Tracker(other.straws(), other.strawProperties(),other.name().c_str()) {
// copy all the G4 variables by hand.  These should never be needed by this copy
    _z0 = other._z0;
    _rOut = other._rOut;
    _envelopeMaterial = other._envelopeMaterial;
    _wallMaterialName = other._wallMaterialName;
    _outerMetalMaterial = other._outerMetalMaterial;
    _innerMetal1Material = other._innerMetal1Material;
    _innerMetal2Material = other._innerMetal2Material;
    _gasMaterialName = other._gasMaterialName;
    _wireMaterialName = other._wireMaterialName;
    _wirePlateMaterial = other._wirePlateMaterial;
    _mother = other._mother;
    _innerTrackerEnvelopeParams = other._innerTrackerEnvelopeParams;
    _planeEnvelopeParams = other._planeEnvelopeParams;
    _panelEnvelopeParams = other._panelEnvelopeParams;
    _supportModel = other._supportModel;
    _supportParams = other._supportParams;
    _panelZOffset = other._panelZOffset;
    _supportStructure = other._supportStructure;
    _manifoldHalfLengths = other._manifoldHalfLengths;
    _envelopeInnerRadius = other._envelopeInnerRadius;
    _panelEB = other._panelEB;
    // other variables: these should move into the functional constructor from straws FIXME!
    _planeExists = other._planeExists;
    _origin = other._origin;
  }

  // G4 accessors: deprecated: FIXME!
  TubsParams Tracker::strawOuterTubsParams(StrawId const& id) const {
    return TubsParams ( 0., strawOuterRadius(), getStraw(id).halfLength() );
  }
  TubsParams Tracker::strawWallMother(StrawId const& id)      const {
    return TubsParams( strawInnerRadius(), 
	strawOuterRadius(), getStraw(id).halfLength() );
  }
  TubsParams Tracker::strawWallOuterMetal(StrawId const& id)  const {
    double rIn = strawOuterRadius() - outerMetalThickness();
    return TubsParams( rIn, strawOuterRadius() , 
		       getStraw(id).halfLength() );
  }
  TubsParams Tracker::strawWallInnerMetal1(StrawId const& id) const {
    double rIn  = strawInnerRadius() + innerMetal2Thickness();
    double rOut = rIn + innerMetal1Thickness();
    return TubsParams (rIn,rOut,
		       getStraw(id).halfLength());
  }
  TubsParams Tracker::strawWallInnerMetal2(StrawId const& id) const {
    double rIn  = strawInnerRadius();
    double rOut = rIn + innerMetal2Thickness();
    return TubsParams (rIn,rOut,
		       getStraw(id).halfLength());
  }
  TubsParams Tracker::strawWireMother(StrawId const& id) const {
    return TubsParams (0.,wireRadius(),
		       getStraw(id).halfLength());
  }
  TubsParams Tracker::strawWirePlate(StrawId const& id) const {
    double rIn = wireRadius() - wirePlateThickness();
    double rOut = wireRadius();
    return TubsParams (rIn,rOut,
		       getStraw(id).halfLength());
  }

} // namespace mu2e
