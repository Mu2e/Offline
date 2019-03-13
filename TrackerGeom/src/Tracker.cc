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
    _stations = other._stations;
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
    _allStraws2 = other._allStraws2;
    _allStraws2_p = other._allStraws2_p;
    _strawExists2 = other._strawExists2;

    // now need to reseat a lot of pointers

    // reseat the non-dense Straw* array
    for(size_t i=0; i<_allStraws2_p.size(); i++ ) {
      _allStraws2_p[i] = &_allStraws2[0] + 
	      (other._allStraws2_p[i] - &other._allStraws2[0]);
    }

    // there are actually two copies of Planes and Panels
    // one held by Stations and one held by Tracker
    // update the Station ones first
    for(size_t istation=0; istation<_stations.size(); istation++) {
      Station station = _stations[istation];
      Station const& ostation = other._stations[istation];

      for(size_t iplane=0; iplane<station._planes.size(); iplane++) {
	Plane& plane = station._planes[iplane];
	plane.fillPointers(this);
	Plane const& oplane = ostation._planes[iplane];

	for(size_t ipanel=0; ipanel<plane._panels.size(); ipanel++) {
	  Panel panel = plane._panels[ipanel];
	  Panel const& opanel = oplane._panels[ipanel];

	  for(size_t i=0; i<panel._straws2_p.size(); i++ ) {
	    panel._straws2_p[i] = &_allStraws2[0] + 
	      (opanel._straws2_p[i] - &other._allStraws2[0]);
	  } // straws

	} // panels

      } // planes

    } // stations

    // Now update the Planes and Panels held by Tracker
    for(size_t iplane=0; iplane<_planes.size(); iplane++) {
      Plane& plane = _planes[iplane];
      plane.fillPointers(this);
      Plane const& oplane = other._planes[iplane];

      for(size_t ipanel=0; ipanel<plane._panels.size(); ipanel++) {
	Panel panel = plane._panels[ipanel];
	Panel const& opanel = oplane._panels[ipanel];

	for(size_t i=0; i<panel._straws2_p.size(); i++ ) {
	  panel._straws2_p[i] = &_allStraws2[0] + 
	    (opanel._straws2_p[i] - &other._allStraws2[0]);
	} // straws

      } // panels

    } // planes

  }

  /*
  Tracker& Tracker::operator=(const Tracker& other) {
    return *this = Tracker(other);
  }
  */

  void Tracker::fillPointers () const{
    for ( size_t i=0; i<StrawId::_nplanes; ++i){
      _planes[i].fillPointers(this);
    }
  }

} // namespace mu2e
