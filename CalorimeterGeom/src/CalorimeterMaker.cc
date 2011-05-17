//
// Make a Calorimeter.
//
// $Id: CalorimeterMaker.cc,v 1.17 2011/05/17 15:35:59 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:35:59 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>

//
// Mu2e includes
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/CalorimeterMaker.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "BeamlineGeom/inc/Beamline.hh"

// Framework include files
#include "cetlib/exception.h"

//
// other includes
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"
using namespace std;




//
// naming conventions:  
// a) a vane is one of the sets of crystals; 
// b) a "z-slice" is one set of crystals in a vane with long side aligned, at constant
//    z in an ideal array;
// c) an "r-slice" is one set of crystals in a vane at the same distance
//    from the z-axis


namespace mu2e{

    CalorimeterMaker::CalorimeterMaker( SimpleConfig const& config)
    {
      _calo = auto_ptr<Calorimeter>(new Calorimeter());
      
      _calo->_nVane      = config.get<int>   ("calorimeter.numberOfVanes");
      _calo->_crystalHW  = config.getDouble("calorimeter.crystalHalfTrans");
      _calo->_crystalHL  = config.getDouble("calorimeter.crystalHalfLong");
      _calo->_nCrystalR  = config.get<int>   ("calorimeter.nCrystalRSlices");   
      _calo->_nCrystalZ  = config.get<int>   ("calorimeter.nCrystalZSlices");   
      _calo->_rMin       = config.getDouble("calorimeter.rInscribed");

      _calo->_wrapperHalfThickness = config.getDouble("calorimeter.crystalWrapperHalfThickness");
      _calo->_roHalfTrans          = config.getDouble("calorimeter.crystalReadoutHalfTrans");
      _calo->_roHalfThickness      = config.getDouble("calorimeter.crystalReadoutHalfThickness");

      GeomHandle<Beamline> beamg;
      double solenoidOffset = beamg->solenoidOffset();
      CLHEP::Hep3Vector center = config.getHep3Vector("calorimeter.calorimeterCenter");
      _calo->_origin = CLHEP::Hep3Vector(-solenoidOffset,0,center.z());

      // Check number of readouts
      int nRO = config.get<int>("calorimeter.crystalReadoutChannelCount");
      if( ! (nRO==1 || nRO==2 || nRO==4) ) {
	throw cet::exception("CaloGeom")
	  << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";
      }
      _calo->_nROPerCrystal = nRO;

      // Check size of readouts
      if( nRO==1 ) {
	if( _calo->_roHalfTrans > _calo->_crystalHW ) {
	  throw cet::exception("CaloGeom")
	    << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHalfTrans.\n";
	}
      } else {
	if( _calo->_roHalfTrans > 0.5*_calo->_crystalHW ) {
	  throw cet::exception("CaloGeom")
	    << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";
	}
      }

      _calo->_nonUniformity = config.getDouble("calorimeter.crystalNonUniformity",0.0);
      _calo->_timeGap       = config.getDouble("calorimeter.timeGap",100.0);
      _calo->_electronEdep  = config.getDouble("calorimeter.electronDepositionAPD",1000.0);
      _calo->_electronEmin  = config.getDouble("calorimeter.electronMinEnergyAPD",0.1);

      // Create vanes
      MakeVanes();

      //
      //make sure above information is consistent
      //CheckIt();

      //
      //do the work
      BuildIt();
    }


    CalorimeterMaker::~CalorimeterMaker() {}



    void CalorimeterMaker::MakeVanes()
    {
      /*
      // Calculate position of DS3 - mother volume of calorimeter

      double rTorus         = beamg->getTS().torusRadius();
      double ts5HalfLength  = beamg->getTS().getTS5().getHalfLength();

      // This probably should go to BeamlineGeom
      double ds2HalfLength     = _config->getDouble("toyDS2.halfLength");
      double ds3HalfLength     = _config->getDouble("toyDS3.halfLength");
      double ds2Z0     = rTorus + 2.*ts5HalfLength + ds2HalfLength;
      double ds3Z0     = ds2Z0  + ds2HalfLength    + ds3HalfLength;

      CLHEP::Hep3Vector posDS3 (-solenoidOffset, 0., ds3Z0);
      */

      // Local reference frame of the vane:
      // Z is along Mu2e z
      // Y is radial towards large radius
      // X is transverse to vane plane

      double dR = _calo->_crystalHW * _calo->_nCrystalR;
      double dZ = _calo->_crystalHW * _calo->_nCrystalZ;
      double radius = _calo->_rMin + dR;
      double dphi = 2*CLHEP::pi/_calo->_nVane;

      for( int i=0; i<_calo->_nVane; ++i ) {

	_calo->_vanes.push_back(Vane(i));
	Vane & v = _calo->_vanes.back();

	v._size = CLHEP::Hep3Vector( _calo->_crystalHL+_calo->_roHalfThickness, dR, dZ );

	double phi = -CLHEP::pi + i*dphi;
	v._origin = CLHEP::Hep3Vector( _calo->_origin.x()+radius*cos(phi),
				       _calo->_origin.y()+radius*sin(phi),
				       _calo->_origin.z() );
	v._originLocal = CLHEP::Hep3Vector(radius*cos(phi),radius*sin(phi),0);
				      
	v._rotation = CLHEP::HepRotation::IDENTITY;
	v._rotation *= CLHEP::HepRotationZ(-CLHEP::pi/2 - i*dphi);

      }

    }

    void CalorimeterMaker::BuildIt()
    {

      // Create material

      //MakeVanes();
      //FillNearestNeighbours();
      //MakeCalorimeter();

      return;
    }

}

