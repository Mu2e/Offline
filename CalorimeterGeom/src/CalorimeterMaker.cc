//
// Make an Calorimeter.
//

//
// C++ includes
#include <iostream>
#include <iomanip>
#include <cmath>

//
// Mu2e includes
#include "Calorimeter/inc/CalorimeterMaker.hh"
#include "Calorimeter/inc/Calorimeter.hh"
#include "Calorimeter/inc/CrystalId.hh"
#include "Calorimeter/inc/Crystal.hh"

// Framework include files
#include "FWCore/Utilities/interface/Exception.h"


//
// other includes
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"


using CLHEP::Hep3Vector;
using CLHEP::HepRotationY;
using CLHEP::HepRotationZ;

using namespace std;

void crystalPrinter( const Crystal& s){
  cout << s.Id() << endl;
}

void crystalPrinter2( const Crystal* s, int& i){
  cout << s->Id() <<  " | " 
       << s->hack << " " 
       << ++i << endl;
}

void layerPrinter( const Layer& l){
  cout << "    Layer: " << l.Id() << endl;
}

void devicePrinter( const Device& d){
  cout << "  Device: " << d.Id() << endl;
}

CalorimeterMaker::CalorimeterMaker(int nVanes,
			     std::vector<LayerInfo> vaneInfo,
			     double r0,
			     double rInscribed,
    			     double caloZLength,
			     double halfSide,
			     Hep3Vector center,
			     double phi0
		   ):
  _nVanes(nVanes),
  _vaneInfo(vaneInfo),
  _r0(r0),
  _rInscribed(rInscribed),
  _caloZLength(caloZLength),
  _halfSide(halfSide),
  _center(center),
  _phi0(phi0)
{
  // Make sure that information is self consisent.
  CheckFit();

  // Make an empty Calorimeter.
  _calo = new Calorimeter();

  // We require that push_back on _allCrystals does not invalidate iterators. 
  // So size it explicitly at the beginning.
  _calo->_allCrystals.reserve(totalCrystals());

  // Set some information.
  _calo->_r0 = r0;
  _calo->_rInscribed = rInscribed;

  // Describe the different types of crystals in the system.
  MakeDetails();

  // Populate all information except the nearest neighbour info.
  MakeVanes();

  // Populate the nearest neighbour info.
  FillNearestNeighbours();

}
  
CalorimeterMaker::~CalorimeterMaker (){}


// Check that the requested number of crystals actually fits.
void CalorimeterMaker::CheckFit (){

  // Check that the layers fit radially
  vector<LayerInfo>::size_type nlayers = _vaneInfo.size();
  double radDiff = _r0-_rInscribed;
  double vaneRad = 2.*_halfSide*nlayers;

  if ( _rInscribed+vaneRad > _r0+radDiff )
    throw cms::Exception("GEOM")
      << "Number of layers with this halfSide do not fit: " << "\n"
      <<"nlayers:      " << nlayers << "\n"
      << "halfSide:     " << _halfSide << " cm\n"
      << "vaneRad: " << vaneRad << " cm is > maxRadius: " 
      << 2.*radDiff << " cm";
  

  // Check that the crystals fit in Z
  for( vector<LayerInfo>::size_type i=0;
       i<_vaneInfo.size(); ++i){
    int nCrystals = _vaneInfo[i]._nCrystals;
    double vaneZ = 2.*_halfSide*nCrystals;
    
    if ( vaneZ > _caloZLength ){
      throw cms::Exception("GEOM")
           << "Number of crystals with this halfSide do not fit: " << "\n"
           << "Layer:        " << i << "\n"
           << "nCrystals:    " << nCrystals << "\n"
           << "halfSide:     " << _halfSide << " cm\n"
	   << "vaneZ: " << vaneZ << " cm is > caloZLength: "
	   << _caloZLength << "cm";
    }
  }
}

// Total number of crystals in the full system.
int CalorimeterMaker::totalCrystals()const{
  int n(0);
  for ( vector<LayerInfo>::size_type i = 0;
	i< _vaneInfo.size(); ++i ){
    n += _vaneInfo[i]._nCrystals;
  }
  n *= _nVanes;
  return n;
}


// More info to be added later.
void CalorimeterMaker::MakeDetails(){
 
  _calo->_crystalDetail.push_back( CrystalDetail( 0, _halfSide*4, _halfSide, _halfSide));
  _calo->_crystalDetail.push_back( CrystalDetail( 0, _halfSide, _halfSide*4, _halfSide));

}

// Build the vanes according to the correct geometry
void CalorimeterMaker::MakeVanes(){

  vector<Crystal>& allCrystals = _calo->_allCrystals;

  vector<LayerInfo>::size_type nlayers = _vaneInfo.size();
  int halfnlayers = nlayers/2;
  
  vector<Device>& dev = _calo->_devices;

  // Loop over the 4 vanes
  for (int vaneNum = 0; vaneNum != 1; ++vaneNum){
    dev.push_back( Device(vaneNum) );
    Device& vane = dev[vaneNum];

    // Initial wire in the X direction
    Hep3Vector baseWire(1.,0.,0.);

    vector<Layer>& layers = vane._layers;
    layers.reserve(nlayers);

    // Define the rotation from the canonical position.
    // 90 degrees between each vane.
    double angle = 2.0*vaneNum*M_PI/_nVanes  + _phi0;
    HepRotationZ RZ(angle);
    
    // Direction for all "wires" in this vane.
    // Imagine each crystal has a wire through the center.
    Hep3Vector wireDir = RZ*baseWire;
    
    DeviceId devId(vaneNum);
    
    // Loop over each layer in the vane
    for ( vector<LayerInfo>::size_type j=0; 
          j<nlayers; ++j ){

      LayerInfo& linfo = _vaneInfo[j];
      int n = linfo._nCrystals;

      // Center of first wire in this layer.
      // Assume the XY origin in in the center of the vane
      double x0 = 0.;
      double y0 = (nlayers%2 == 1)? 
        (-1.*double(j)+halfnlayers) :
        (-1.*double(j)+halfnlayers-0.5);
      
      // Set origin of this wire and offset to the next in the layer
      Hep3Vector origin = Hep3Vector( x0, (y0*2.*_halfSide)-_r0, 0.);
      Hep3Vector delta  = Hep3Vector( 0., 0., 2*_halfSide);
      
      // Store the information
      LayerId lid(devId,j);
      layers.push_back( Layer(lid, linfo._nCrystals, origin, delta));
      
      vector<CrystalIndex>& indices  = layers[j]._indices;
      indices.reserve(n);

      vector<const Crystal*>& crystals = layers[j]._crystals;

      // Choose which CrystalDetail to save according to vaneNum
      // See void CalorimeterMaker::MakeDetails above
      int detailIndex = 0;
      if (devId == 1 || devId == 3){
        detailIndex = 1;
      }
      CrystalDetail* detail = &_calo->_crystalDetail[detailIndex];
      
      // Loop over each crystal
      for ( int m=0; m<linfo._nCrystals; ++m ){

	// Set the center of this crystal and rotate it properly.
        Hep3Vector p = origin + m*delta;  
        Hep3Vector q = RZ*p;

        CrystalIndex index = allCrystals.size();

        // This is the operation that required the call to .reserve in the c'tor.
        allCrystals.push_back( Crystal( CrystalId(lid,m), 
      			    index, q, detail, detailIndex, wireDir) );

        // Fill the other containers.
        indices.push_back( index );	
        crystals.push_back( &allCrystals[index] );
      }      
    }
  }
}

// Check that each layer has the same number of crystals
void CalorimeterMaker::CheckVaneConsistency(){
  if ( _vaneInfo.size()<2 ){
    cerr << "Vanes with fewer than 2 layers are not supported." << endl;
    exit(-1);
    // throw;
  }
  for ( vector<LayerInfo>::size_type i = 1;
	i<_vaneInfo.size(); ++i ){
    unsigned int j = i-1;
    int delta = _vaneInfo[i]._nCrystals - _vaneInfo[j]._nCrystals;
    if (delta != 0 ){
      cerr << "This version only supports vanes in which crystals/layer does not change."
	   << endl;
      exit(-1);
      // throw here.
    }
  }

}


//
// Calculate a crystal's nearest neighbors in 2 passes.
// Pass 1: Fill the _nearestById container.
// Pass 2: Fill the _nearest and _nearestByIndex containers.
//
void CalorimeterMaker::FillNearestNeighbours(){
  
  // Build the nearest neighbour info for each crystal.
  for ( vector<Crystal>::iterator i=_calo->_allCrystals.begin(), 
	  e=_calo->_allCrystals.end();
	i!=e; 
	++i){
    
    // For readability.
    Crystal& cryst = *i;

    // Get references to the device and layer that hold this crystal.
    const DeviceId& devid = cryst.Id().getDeviceId();
    const LayerId& layid  = cryst.Id().getLayerId();
    const Device& device  = _calo->getDevice(devid);
    const Layer& layer    = _calo->getLayer(layid);

    int jL = cryst.Id().getLayer();
    int jC = cryst.Id().getCrystal();

    // Number of crystals in this layer.
    int ncrystals = layer.nCrystals();

    // Number of layers in this device.
    int nlayers = device.nLayers();

    // Add crystals in the same layer.  Deal with edge cases.
    if ( jC == 0 ){
      cryst._nearestById.push_back(CrystalId(devid,jL,jC+1));
    } else if ( jC == ncrystals-1 ){
      cryst._nearestById.push_back(CrystalId(devid,jL,jC-1));
    } else {
      cryst._nearestById.push_back(CrystalId(devid,jL,jC-1));
      cryst._nearestById.push_back(CrystalId(devid,jL,jC+1));
    }
      
    // Add crystals one layer inward one layer, unless already in innermost layer. 
    // Deal with edge cases.
    if ( jL != 0 ){
      if ( jC == 0 ){
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC));
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC+1));
	} else if ( jC == ncrystals-1 ) {
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC));
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC-1));
	}else{
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC));
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC-1));
	cryst._nearestById.push_back(CrystalId(devid,jL-1,jC+1));
	}
      }

    // Add crystals from one layer outward, unless already in outermost layer. 
    // Deal with edge cases.
    if ( jL != nlayers-1){
      if ( jC == 0 ){
	cryst._nearestById.push_back(CrystalId(devid,jL+1,jC));
        cryst._nearestById.push_back(CrystalId(devid,jL+1,jC+1));
      } else if ( jC == ncrystals-1 ){
        cryst._nearestById.push_back(CrystalId(devid,jL+1,jC));
	cryst._nearestById.push_back(CrystalId(devid,jL+1,jC-1));
      } else{
	cryst._nearestById.push_back(CrystalId(devid,jL+1,jC));
        cryst._nearestById.push_back(CrystalId(devid,jL+1,jC+1));
	cryst._nearestById.push_back(CrystalId(devid,jL+1,jC-1));
      }
    }
    
  }

  // Second pass.  Fill the other two containers.
  
  // Fill nearest neighbour indices and pointers from the NN Ids.
  for ( vector<Crystal>::iterator i= _calo->_allCrystals.begin(), 
	  e= _calo->_allCrystals.end();
	i!=e; 
	++i){
    vector<const Crystal *>& byPtr = i->_nearest;
    vector<CrystalId>& byId        = i->_nearestById;
    vector<CrystalIndex>& byIndex  = i->_nearestByIndex;
    
    byPtr.clear();
    byIndex.clear();
    
    for ( vector<CrystalId>::iterator j=byId.begin(), je=byId.end();
	   j != je; ++j){
      const CrystalId& id = *j;
      const Crystal& crystal = _calo->getCrystal(id);
      byPtr.push_back( &crystal);
      byIndex.push_back( crystal.Index() );
    }
  }

}

