//
// Construct and return an LTracker.
//
//
// $Id: LTrackerMaker.cc,v 1.5 2010/02/06 22:20:28 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/02/06 22:20:28 $
//
// Original author Rob Kutschke
//

#include <iostream>
#include <iomanip>
#include <cmath>


// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "LTrackerGeom/inc/LTrackerMaker.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/StrawId.hh"
#include "LTrackerGeom/inc/Straw.hh"
#include "Mu2eUtilities/inc/for_all.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"

#ifndef __CINT__ 

using namespace CLHEP;

using namespace std;

namespace mu2e {

  void strawPrinter( const Straw& s){
    cout << s.Id() << endl;
  }

  void strawPrinter2( const Straw* s, int& i){
    cout << s->Id() <<  " | " 
	 << s->hack << " " 
	 << ++i << endl;
  }

  void strawHacker( Straw* s, int& i){
    s->hack = 2;
  }
  
  void layerPrinter( const Layer& l){
    cout << "    Layer: " << l.Id() << endl;
  }
  
  void sectorPrinter( const Sector& s){
    cout << "  Sector: " << s.Id() << endl;
  }
  
  void devicePrinter( const Device& d){
    cout << "  Device: " << d.Id() << endl;
  }
  
  
  LTrackerMaker::LTrackerMaker(int nSides,
			       std::vector<LayerInfo> sideInfo,
			       std::vector<LayerInfo> vaneInfo,
			       double r0,
			       double halfLength,
			       double radius,
			       Hep3Vector center,
			       double phi0,
			       double tilt,
			       Hep3Vector vaneOffset
			       ):
    _nSides(nSides),
    _sideInfo(sideInfo),
    _vaneInfo(vaneInfo),
    _r0(r0),
    _halfLength(halfLength),
    _strawRadius(radius),
    _center(center),
    _phi0(phi0),
    _tilt(tilt),
    _vaneOffset(vaneOffset)
  {
    // Do the work.
    BuildIt();
  }

  // Constructor that gets information from the config file instead of
  // from arguments.
  LTrackerMaker::LTrackerMaker( SimpleConfig const& config):
    _center(),
    _sideInfo(),
    _vaneInfo(){

    _nSides       = config.getInt("ltracker.nSides");
    _r0           = config.getDouble("ltracker.r0");
    _z0           = config.getDouble("ltracker.z0");
    _halfLength   = config.getDouble("ltracker.zHalfLength");
    _rOut         = config.getDouble("ltracker.rOut");
    _strawRadius  = config.getDouble("ltracker.rStrawOut");
    _phi0         = config.getDouble("ltracker.phi0");
    _tilt         = config.getDouble("ltracker.tilt");
    _strawThick   = config.getDouble("ltracker.strawThickness");
    _rwire        = config.getDouble("ltracker.rWire");
    _carbonThick  = config.getDouble("ltracker.carbonThick");
    _vaneOffset   = config.getHep3Vector("ltracker.vaneOffset");
    _fillMaterial = config.getString("ltracker.fillMaterial");

    config.getVectorString("ltracker.strawMaterials0", _strawMaterialNames0, 3);
    config.getVectorString("ltracker.strawMaterials1", _strawMaterialNames1, 3);
    
    vector<int> nStrawsSide, nStrawsVane;
    config.getVectorInt("ltracker.nStrawsSide", nStrawsSide);
    config.getVectorInt("ltracker.nStrawsVane", nStrawsVane);
    
    LayerInfo::Stype edge = LayerInfo::nonconductive;
    LayerInfo::Stype mid  = LayerInfo::conductive;
    
    for ( int i=0; i<nStrawsSide.size(); ++i ){
      LayerInfo::Stype stype = ( i == 0 || i == nStrawsSide.size()-1 ) ? edge : mid;
      _sideInfo.push_back( LayerInfo( nStrawsSide[i], stype) );
    }
    
    for ( int i=0; i<nStrawsVane.size(); ++i ){
      LayerInfo::Stype stype = ( i == 0 || i == nStrawsVane.size()-1 ) ? edge : mid;
      _vaneInfo.push_back( LayerInfo( nStrawsVane[i], stype) );
    }

    // Do the real work.
    BuildIt( );
  }


  
  LTrackerMaker::~LTrackerMaker (){}

  void LTrackerMaker::BuildIt(){

    _phiHalf = M_PI/_nSides;
    _tphiHalf = tan(_phiHalf);
    _cphiHalf = cos(_phiHalf);
    _sphiHalf = sin(_phiHalf);
    
    CheckFit();

    _ltt = auto_ptr<LTracker>(new LTracker());
    
    _ltt->_r0         = _r0;
    _ltt->_z0         = _z0;
    _ltt->_rOut       = _rOut;
    _ltt->_rInscribed = _r0-_strawRadius*(1.+sqrt(3.));
    
    _ltt->_fillMaterial = _fillMaterial;

    MakeDetails();
    MakeSides();
    MakeVanes();

    // Fill the pointers inside the Layer level.
    _ltt->FillPointers1();

    /*
      // Testing.
      cout << "Straw Printer: " << endl;
      _ltt->forAllStraws(strawPrinter);
      
      cout << "Layer Printer: " << endl;
      _ltt->forAllLayers(layerPrinter);
      
      cout << "Sector Printer: " << endl;
      _ltt->forAllSectors(sectorPrinter);
      
      cout << "Device Printer: " << endl;
      _ltt->forAllDevices(devicePrinter);
    */
    //const Layer& l0 = _ltt->getLayer( LayerId(LTracker::wedge,0,0) );
    //const vector<const Straw*>& v = l0.getStraws();
    
    //int nn(-1);
    //  for_all ( v, strawHacker, nn);
    //  for_all ( v, strawPrinter2, nn);

    //int mm(-1);
    //for_all ( v, test, mm);


    // Must do this after the LTracker has been built since this 
    // computes pointers to things that must remain in  static
    // memory locations.
    
    FillNearestNeighbours();
    _ltt->FillPointers2();

  }

  // Check that the requested number of straws actually fits.
  void LTrackerMaker::CheckFit (){

    for( vector<LayerInfo>::size_type i=0;
	 i<_sideInfo.size(); ++i){
      LayerInfo& lay = _sideInfo[i];
      
      double f = 1. + double(i)*sqrt(3.);

      double r = _r0-_strawRadius*(1.+sqrt(3.));
      
      double maxRadius = (r*_tphiHalf)/
	( lay._nStraws - 1. + _cphiHalf - (f-_sphiHalf)*_tphiHalf);
      
      if ( _strawRadius > maxRadius ){
	throw cms::Exception("GEOM")
	  << "Straw radius too big for straws to fit in the defined space: \n"
	  << "Layer:      " << i << "\n"
	  << "rInscribed: " << r << "\n"
	  << "nsides:     " << _nSides << "\n"
	  << "nstraws:    " << lay._nStraws << "\n"
	  << "radius:     " << _strawRadius << "\n"
	  << "maxradius:  " << maxRadius
	<< "\n";
      }
      
    }
  }


  int LTrackerMaker::totalStraws()const{
    int n(0);
    for ( vector<LayerInfo>::size_type i = 0;
	  i< _sideInfo.size(); ++i ){
      n += _sideInfo[i]._nStraws;
    }
    n *= _nSides;
    return n;
  }
  
  void LTrackerMaker::MakeDetails(){
        
    _ltt->_strawDetail.push_back( 
	 StrawDetail( _strawMaterialNames0, _strawRadius, _strawThick, _halfLength, _rwire) );
    _ltt->_strawDetail.push_back( 
         StrawDetail( _strawMaterialNames1, _strawRadius, _strawThick, _halfLength, _rwire) );
  }

  void LTrackerMaker::MakeSides(){

    CheckSideConsistency();

    deque<Straw>& allStraws = _ltt->_allStraws;

    // Needed if we change to a vector.
    //  allStraws.reserve(totalStraws());

    vector<LayerInfo>::size_type nlayers = _sideInfo.size();
    int halfnlayers = nlayers/2;

    vector<Device>& dev = _ltt->_devices;

    dev.push_back( Device(LTracker::wedge) );
    Device& sides = dev[LTracker::wedge];
    vector<Sector>& sectors = sides._sectors;
    sectors.reserve(_nSides);

    const double root3(sqrt(3.0));
    double yoffset(root3*_strawRadius);

    // Rotation in the plane of the sector.
    HepRotationY RY(-_tilt);

    Hep3Vector baseWire(0.,0.,1.);

    for ( int isec=0; isec<_nSides; ++isec ){

      sectors.push_back(Sector(SectorId(LTracker::wedge,isec)));
      Sector& sec = sectors.back();
      
      vector<Layer>& layers = sec._layers;
      layers.reserve(nlayers);

      // Define the rotation from the canonical position.
      double angle = isec*2.*M_PI/_nSides - 3*M_PI/_nSides + _phi0;
      HepRotationZ RZ(angle);

      // Offset to place the center of the sector.
      double angle2 = isec*2.*M_PI/_nSides + M_PI/_nSides + _phi0;
      Hep3Vector sectorOffset(_r0*cos(angle2),_r0*sin(angle2),0.);

      // Wire direction for all wires in this sector.
      Hep3Vector wireDir = RZ*(RY*baseWire);
    
      SectorId secId(LTracker::wedge,isec);

      int maxStrawsThisSector(0);

      for ( vector<LayerInfo>::size_type j=0; 
	    j<nlayers; ++j ){
	
	LayerInfo& linfo = _sideInfo[j];
	int n = linfo._nStraws;

	maxStrawsThisSector = ( n > maxStrawsThisSector ) ? n : maxStrawsThisSector;
	
	// For a sector centered at the origin, compute location of wire 0 
	// and offset between wires within this layer.
	double x0 = (n-1)*_strawRadius;
	double y0 = ( nlayers%2 == 1)? 
	  (double(j)-halfnlayers)*yoffset :
	  (double(j)-halfnlayers+0.5)*yoffset;
	
	Hep3Vector origin = Hep3Vector( x0,   y0, 0.);
	Hep3Vector delta  = Hep3Vector( -2.*_strawRadius, 0., 0.);
      
	LayerId lid(secId,j);
	layers.push_back( Layer(lid, n, origin, delta));
	
	//vector<const Straw*>& straws = layers[j]._straws;
	vector<StrawIndex>& indices  = layers[j]._indices;

	int detailIndex = linfo._strawType;
	StrawDetail* detail = &_ltt->_strawDetail[detailIndex];

	// Save the base position of straw 0 in the sector.
	sec._basePosition.push_back(origin);

	for ( int is=0; is<linfo._nStraws; ++is ){

	  // Position in sector centered at the origin.
	  Hep3Vector p = origin + is*delta;

	  // Index into master container.
	  StrawIndex index = StrawIndex(allStraws.size());
	  
	  // Final position with rotation and translation.
	  Hep3Vector q = RZ*(RY*p) + sectorOffset;

	  // This operation must not invalidate pointers to
	  // previous elements in allStraws.
	  allStraws.push_back( Straw( StrawId(lid,is), 
	        index, q, detail, detailIndex, wireDir) );
	  indices.push_back( index );
	  
	  //This is why two lines back must not invalidate pointers.
	  //straws.push_back( &allStraws[index] );

	} // end loop over straws

      } // end loop over layers

      // Offset between wires in base position.
      sec._baseDelta = Hep3Vector( -2.*_strawRadius, 0., 0.);

      // Scale the box to be slightly larger than it needs to be
      const double fac = 1.0001;

      // Half-dimensions of a box that holds this layer.
      sec._boxHalfLengths.push_back( fac * (maxStrawsThisSector   * _strawRadius) );
      sec._boxHalfLengths.push_back( fac * (1.+0.5*(nlayers-1)*root3)* _strawRadius );
      sec._boxHalfLengths.push_back( fac * _halfLength );

      // Descriptions of the rotations and placement of the sector.
      sec._boxRyAngle = -_tilt;
      sec._boxRzAngle = angle;
      sec._boxOffset  = sectorOffset;


    } // end loop over sectors

  } // end make sides.


void LTrackerMaker::MakeVanes(){
  
  deque<Straw>& allStraws = _ltt->_allStraws;

  vector<LayerInfo>::size_type nlayers = _sideInfo.size();
  int halfnlayers = nlayers/2;
  
  vector<Device>& dev = _ltt->_devices;

  dev.push_back( Device(LTracker::vane) );
  Device& vanes = dev[LTracker::vane];
  vector<Sector>& sectors = vanes._sectors;
  sectors.reserve(_nSides);

  const double root3(sqrt(3.0));
  double yoffset(root3*_strawRadius);

  HepRotationY RY(-_tilt);

  Hep3Vector baseWire(0.,0.,1.);

  for ( int isec=0; isec<_nSides; ++isec ){

    sectors.push_back(Sector(SectorId(LTracker::vane,isec)));
    Sector& sec = sectors.back();
   
    vector<Layer>& layers = sec._layers;
    layers.reserve(nlayers);

    // Define the rotation from the canonical position.
    double angle = 2.*isec*M_PI/_nSides  + _phi0;
    HepRotationZ RZ(angle);
    
    // Wire direction for all wires in this sector.
    Hep3Vector wireDir = RZ*(RY*baseWire);
    
    SectorId secId(LTracker::vane,isec);

    int maxStrawsThisSector(0);
    
    for ( vector<LayerInfo>::size_type j=0; 
	  j<nlayers; ++j ){

      LayerInfo& linfo = _vaneInfo[j];
      
      int n = linfo._nStraws;

      maxStrawsThisSector = ( n > maxStrawsThisSector ) ? n : maxStrawsThisSector;
      
      // Center of first wire in this layer.
      double x0 = -(n-1)*_strawRadius;
      double y0 = ( nlayers%2 == 1)? 
	(double(j)-halfnlayers)*yoffset :
	(double(j)-halfnlayers+0.5)*yoffset;

      // For testing purposes
      //y0 += 4*isec + 2;
      
      // Rotate origin and offset to this sector.
      Hep3Vector origin = Hep3Vector( x0,   y0, 0.);
      Hep3Vector delta  = Hep3Vector( 2.*_strawRadius, 0., 0.);
      
      LayerId lid(secId,j);
      layers.push_back( Layer(lid, linfo._nStraws, origin, delta));
      
      vector<StrawIndex>& indices  = layers[j]._indices;
      indices.reserve(n);
      
      int detailIndex = linfo._strawType;
      StrawDetail* detail = &_ltt->_strawDetail[detailIndex];

      // Save the base position of straw 0 in the sector.
      sec._basePosition.push_back(origin);

      for ( int is=0; is<linfo._nStraws; ++is ){
	Hep3Vector p = origin + is*delta;
	StrawIndex index = StrawIndex(allStraws.size());

	Hep3Vector q = RZ*(RY*p + _vaneOffset);

	// This operation must not invalidate pointers to
	// previous elements in allStraws.
	allStraws.push_back( Straw( StrawId(lid,is), 
				    index, q, detail, detailIndex, wireDir) );
	indices.push_back( index );
	
	// This is why two lines back must not invalidate pointers.
	//straws.push_back( &allStraws[index] );
      } // end loop over straws
      
    } // end loop over layers

    // Offset for constructing straw positions.
    sec._baseDelta = Hep3Vector( 2.*_strawRadius, 0., 0.);

    // Half-dimensions of a box that holds this layer.
    sec._boxHalfLengths.push_back( maxStrawsThisSector   * _strawRadius );
    sec._boxHalfLengths.push_back( (1.+0.5*(nlayers-1)*root3)* _strawRadius );
    sec._boxHalfLengths.push_back( _halfLength );

    // Descriptions of the rotations and placement of the sector.
    sec._boxRyAngle = -_tilt;
    sec._boxRzAngle = angle;
    sec._boxOffset  = RZ*_vaneOffset;

  } // end loop over sectors

}

// Two parts of the code assume that each layer has one more straw than the
// next innermost.
//  1) The nearest neighbour algorithm.
//  2) The code to compute the bounding volumes.
void LTrackerMaker::CheckSideConsistency(){
  if ( _sideInfo.size()<2 ){
    throw cms::Exception("GEOM")
      << "Side sectors with fewer than 2 layers are not supported.\n";
  }
  for ( vector<LayerInfo>::size_type i = 1;
	i<_sideInfo.size(); ++i ){
    unsigned int j = i-1;
    if ( _sideInfo[i]._nStraws != _sideInfo[j]._nStraws + 1){
      throw cms::Exception("GEOM")
	<< "This version only supports sectors in which straws/layer increases by 1.\n";
    }
  }
}


// Two parts of the code assume that straws per layer differ by +/1.
//  1) The nearest neighbour algorithm.
//  2) The code to compute the bounding volumes.
void LTrackerMaker::CheckVaneConsistency(){
  if ( _vaneInfo.size()<2 ){
    throw cms::Exception("GEOM")
      << "Vanes with fewer than 2 layers are not supported.\n";
  }
  for ( vector<LayerInfo>::size_type i = 1;
	i<_vaneInfo.size(); ++i ){
    unsigned int j = i-1;
    int delta = _vaneInfo[i]._nStraws - _vaneInfo[j]._nStraws;
    if ( abs(delta) != 1 ){
      throw cms::Exception("GEOM")
	<< "This version only supports vanes in which straws/layer differs by +/1.\n";
    }
  }

}


//
// This algorithm assumes that the number of straws in adjacent
// layers differs by at most 1.
// 
void LTrackerMaker::FillNearestNeighbours(){

  // Build the nearest neighbour info for each straw.
  for ( deque<Straw>::iterator i=_ltt->_allStraws.begin(), 
	  e=_ltt->_allStraws.end();
	i!=e; 
	++i){
    
    // For readability.
    Straw& str = *i;

    // Get references to the sector and layer that hold this straw.
    const SectorId& secid = str.Id().getSectorId();
    const LayerId& layid  = str.Id().getLayerId();
    const Sector& sector  = _ltt->getSector(secid);
    const Layer& layer    = _ltt->getLayer(layid);

    int jl0 = str.Id().getLayer();
    int js0 = str.Id().getStraw();

    // Number of straws in this layer.
    int nstraws = layer.nStraws();

    // Number of layers in this sector.
    int nlayers = sector.nLayers();

    // Add straws in the same layer.  Deal with edge cases.
    if ( js0 == 0 ){
      str._nearestById.push_back(StrawId(secid,jl0,js0+1));
    } else if ( js0 == nstraws-1 ){
      str._nearestById.push_back(StrawId(secid,jl0,js0-1));
    } else {
      str._nearestById.push_back(StrawId(secid,jl0,js0-1));
      str._nearestById.push_back(StrawId(secid,jl0,js0+1));
    }
      
    // Add straws one layer inward one layer, unless already in innermost layer. 
    // Deal with edge cases.
    if ( jl0 != 0 ){
      int nstrawsInside = sector.getLayer(jl0-1).nStraws();
      if ( js0 == 0 ){
	str._nearestById.push_back(StrawId(secid,jl0-1,js0));
	if ( nstraws<nstrawsInside){
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0+1));
	}
      } else if ( js0 == nstraws-1 ) {
	if ( nstraws<nstrawsInside){
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0));
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0+1));
	}else{
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0-1));
	}
      } else{
	str._nearestById.push_back(StrawId(secid,jl0-1,js0));
	if ( nstraws<nstrawsInside){
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0+1));
	}else{
	  str._nearestById.push_back(StrawId(secid,jl0-1,js0-1));
	}
      }
    }

    // Add straws from one layer outward, unless already in outermost layer. 
    // Deal with edge cases.
    if ( jl0 != nlayers-1){
      int nStrawsOutside = sector.getLayer(jl0+1).nStraws();
      if ( js0 == 0 ){
	str._nearestById.push_back(StrawId(secid,jl0+1,js0));
	if ( nstraws<nStrawsOutside){
	    str._nearestById.push_back(StrawId(secid,jl0+1,js0+1));
	}
      } else if ( js0 == nstraws-1 ){
	if ( nstraws < nStrawsOutside ){
	  str._nearestById.push_back(StrawId(secid,jl0+1,js0));
	  str._nearestById.push_back(StrawId(secid,jl0+1,js0+1));
	}else{
	  str._nearestById.push_back(StrawId(secid,jl0+1,js0-1));
	  }
      } else{
	str._nearestById.push_back(StrawId(secid,jl0+1,js0));
	if ( nstraws < nStrawsOutside ){
	  str._nearestById.push_back(StrawId(secid,jl0+1,js0+1));
	}else{
	  str._nearestById.push_back(StrawId(secid,jl0+1,js0-1));
	}
      }
    }
    
    
  }
}


void LTrackerMaker::FillPointersAndIndices(){

  // Fill the main link
  for ( vector<Device>::iterator idev = _ltt->_devices.begin(),
	  edev = _ltt->_devices.end(); 
	idev != edev;  ++idev ){
    for ( vector<Sector>::iterator isec = idev->_sectors.begin(),
	    esec = idev->_sectors.end();
	  isec != esec; ++isec ){
      for ( vector<Layer>::iterator ilay = isec->_layers.begin(),
	      elay = isec->_layers.end();
	    ilay != elay; ++ilay ){
	ilay->_straws.clear();

	for ( vector<StrawIndex>::iterator istr = ilay->_indices.begin(),
		estr = ilay->_indices.end();
	      istr != estr ; ++istr ){
	  const Straw& str = _ltt->_allStraws[(*istr).asInt()];
	  ilay->_straws.push_back( &str );
	}
      }
    }
  }


}
void LTrackerMaker::FillPointersAndIndices2(){

  // Fill nearest neighbour indices and pointers from the NN Ids.
  for ( deque<Straw>::iterator i=_ltt->_allStraws.begin(), 
	  e=_ltt->_allStraws.end();
	i!=e; 
	++i){
    //    vector<const Straw *>& byPtr= i->_nearest;
    vector<StrawId>& byId = i->_nearestById;
    vector<StrawIndex>& byIndex = i->_nearestByIndex;

    //   byPtr.clear();
    byIndex.clear();

    for ( vector<StrawId>::iterator j=byId.begin(), je=byId.end();
	   j != je; ++j){
      const StrawId& id = *j;
      const Straw& straw = _ltt->getStraw(id);
      //      byPtr.push_back( &straw);
      byIndex.push_back( straw.Index() );
    }
  }

}

} // namespace mu2e

#endif
