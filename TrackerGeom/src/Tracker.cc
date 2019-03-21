//
// Geometry and identifier info about an Tracker.
//
//
// $Id: Tracker.cc,v 1.8 2013/01/07 04:01:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/07 04:01:16 $
//
// Original author Rob Kutschke
//

#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {

  Tracker::Tracker(const Tracker& other) {
    // shallow copy
    _name = other._name;
    _z0 = other._z0;
    _rOut = other._rOut;
    _envelopeMaterial = other._envelopeMaterial;
    _strawInnerRadius = other._strawInnerRadius;
    _strawOuterRadius = other._strawOuterRadius;
    _strawWallThickness = other._strawWallThickness;
    _outerMetalThickness = other._outerMetalThickness;
    _innerMetal1Thickness = other._innerMetal1Thickness;
    _innerMetal2Thickness = other._innerMetal2Thickness;
    _wireRadius = other._wireRadius;
    _wirePlateThickness = other._wirePlateThickness;
    _wallMaterialName = other._wallMaterialName;
    _outerMetalMaterial = other._outerMetalMaterial;
    _innerMetal1Material = other._innerMetal1Material;
    _innerMetal2Material = other._innerMetal2Material;
    _gasMaterialName = other._gasMaterialName;
    _wireMaterialName = other._wireMaterialName;
    _wirePlateMaterial = other._wirePlateMaterial;
    _strawHalfLengths = other._strawHalfLengths;
    _strawActiveHalfLengths = other._strawActiveHalfLengths;
    _allManifolds = other._allManifolds;
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
    _planes = other._planes;
    _panels = other._panels;
    _allStraws = other._allStraws;
    _allStraws_p = other._allStraws_p;
    _strawExists2 = other._strawExists2;

    // now need to complete deep copy by
    // reseating a lot of pointers

    // update the Plane's pointers to its Panels
    // and the Panel's pointers to its Straws
    for(size_t iplane=0; iplane<_planes.size(); iplane++) {
      Plane& plane = _planes[iplane];
      plane.fillPointers(this);
      Plane const& oplane = other._planes[iplane];

      for(size_t ipanel=0; ipanel<plane._panels.size(); ipanel++) {
	Panel panel = *( plane._panels[ipanel] );
	Panel const& opanel = *( oplane._panels[ipanel] );

	plane._panels[ipanel] = &_panels[0] + 
	    (oplane._panels[ipanel] - &other._panels[0]);

	for(size_t i=0; i<panel._straws2_p.size(); i++ ) {
	  panel._straws2_p[i] = &_allStraws[0] + 
	    (opanel._straws2_p[i] - &other._allStraws[0]);
	} // straws

      } // panels

    } // planes


    // reseat the non-dense Straw* array
    for(size_t i=0; i<_allStraws_p.size(); i++ ) {
      _allStraws_p[i] = &_allStraws[0] + 
	      (other._allStraws_p[i] - &other._allStraws[0]);
    }


  }

  void Tracker::fillPointers () const{
    for ( size_t i=0; i<StrawId::_nplanes; ++i){
      _planes[i].fillPointers(this);
    }
  }

  TubsParams Tracker::strawOuterTubsParams(StrawId const& id) const {
    return TubsParams ( 0., strawOuterRadius(), getStrawHalfLength(id.straw()) );
  }
  TubsParams Tracker::strawWallMother(StrawId const& id)      const {
    return TubsParams( strawInnerRadius(), 
		       strawOuterRadius(), getStrawHalfLength(id.straw()) );
  }
  TubsParams Tracker::strawWallOuterMetal(StrawId const& id)  const {
    double rIn = strawOuterRadius() - outerMetalThickness();
    return TubsParams( rIn, strawOuterRadius() , 
		       getStrawHalfLength(id.straw()) );
  }
  TubsParams Tracker::strawWallInnerMetal1(StrawId const& id) const {
    double rIn  = strawInnerRadius() + innerMetal2Thickness();
    double rOut = rIn + innerMetal1Thickness();
    return TubsParams (rIn,rOut,
		       getStrawHalfLength(id.straw()));
  }
  TubsParams Tracker::strawWallInnerMetal2(StrawId const& id) const {
    double rIn  = strawInnerRadius();
    double rOut = rIn + innerMetal2Thickness();
    return TubsParams (rIn,rOut,
		       getStrawHalfLength(id.straw()));
  }
  TubsParams Tracker::strawWireMother(StrawId const& id) const {
    return TubsParams (0.,wireRadius(),
		       getStrawHalfLength(id.straw()));
  }
  TubsParams Tracker::strawWirePlate(StrawId const& id) const {
    double rIn = wireRadius() - wirePlateThickness();
    double rOut = wireRadius();
    return TubsParams (rIn,rOut,
		       getStrawHalfLength(id.straw()));
  }

} // namespace mu2e
