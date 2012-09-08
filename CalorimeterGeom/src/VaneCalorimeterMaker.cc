//
// Make a Vane Calorimeter.
//
// $Id: VaneCalorimeterMaker.cc,v 1.1 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>

//
// Mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"

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

  VaneCalorimeterMaker::VaneCalorimeterMaker( SimpleConfig const& config, double solenoidOffset)
  {
      
      _calo = auto_ptr<VaneCalorimeter>(new VaneCalorimeter());

      _calo->_nVane                = config.getInt   ("calorimeter.numberOfVanes");
      _calo->_crystalHW            = config.getDouble("calorimeter.crystalHalfTrans");
      _calo->_crystalHL            = config.getDouble("calorimeter.crystalHalfLong");
      _calo->_nCrystalR            = config.getInt   ("calorimeter.nCrystalRSlices");
      _calo->_nCrystalZ            = config.getInt   ("calorimeter.nCrystalZSlices");
      _calo->_rMin                 = config.getDouble("calorimeter.rInscribed");

      _calo->_wrapperThickness     = config.getDouble("calorimeter.crystalWrapperThickness",0.0);
      _calo->_shellThickness       = config.getDouble("calorimeter.crystalShellThickness",0.0);

      _calo->_nROPerCrystal        = config.getInt("calorimeter.crystalReadoutChannelCount");
      _calo->_roHalfTrans          = config.getDouble("calorimeter.crystalReadoutHalfTrans");
      _calo->_roHalfThickness      = config.getDouble("calorimeter.crystalReadoutHalfThickness");

      _calo->_nonUniformity        = config.getDouble("calorimeter.crystalNonUniformity",0.0);
      _calo->_timeGap              = config.getDouble("calorimeter.timeGap",100.0);
      _calo->_electronEdep         = config.getDouble("calorimeter.electronDepositionAPD",1000.0);
      _calo->_electronEmin         = config.getDouble("calorimeter.electronMinEnergyAPD",0.1);

      _calo->_apdMeanNoise         = config.getDouble("calorimeter.meanNoiseAPD", 0.0);
      _calo->_apdSigmaNoise        = config.getDouble("calorimeter.sigmaNoiseAPD", 0.03);
      _calo->_lysoLightYield       = config.getDouble("calorimeter.lysoLightYield", 2000.0);
      _calo->_apdQuantumEff        = config.getDouble("calorimeter.quantumEffAPD", 0.68);
      _calo->_lightCollectEffAPD   = config.getDouble("calorimeter.lightCollectEffAPD", 0.11);

      
      CLHEP::Hep3Vector center = config.getHep3Vector("calorimeter.calorimeterCenter");
      _calo->_origin = CLHEP::Hep3Vector(-solenoidOffset,0,center.z());

      
      //make sure above information is consistent
      CheckIt();
            
      // Create vanes
      MakeVanes();

    }


    VaneCalorimeterMaker::~VaneCalorimeterMaker() {}



    void VaneCalorimeterMaker::MakeVanes()
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

      double dX     =  _calo->_crystalHL + _calo->_wrapperThickness + _calo->_roHalfThickness;
      double dR     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR;
      double dZ     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalZ;
      double radius = _calo->_rMin + dR;
      double dphi   = 2*CLHEP::pi/_calo->_nVane;


      for( int i=0; i<_calo->_nVane; ++i ) {

        _calo->_vanes.push_back(Vane(i));
        Vane & v = _calo->_vanes.back();

        v._size = CLHEP::Hep3Vector( dX, dR, dZ );

        double phi = -CLHEP::pi + i*dphi;
        v._origin = CLHEP::Hep3Vector( _calo->_origin.x()+radius*cos(phi),
                                       _calo->_origin.y()+radius*sin(phi),
                                       _calo->_origin.z() );
        v._originLocal = CLHEP::Hep3Vector(radius*cos(phi),radius*sin(phi),0);

        v._rotation = CLHEP::HepRotation::IDENTITY;
        v._rotation *= CLHEP::HepRotationZ(-CLHEP::pi/2 - i*dphi);

      }

    }



  void VaneCalorimeterMaker::CheckIt(void)
  {
      // Check number of readouts
      int nRO = _calo->_nROPerCrystal;
      
      if( ! (nRO==1 || nRO==2 || nRO==4) ) 
        {throw cet::exception("CaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}
      

      // Check size of readouts
      if( nRO==1 ) {
         if( _calo->_roHalfTrans > _calo->_crystalHW ) 
           {throw cet::exception("CaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHalfTrans.\n";}
        
      } else {
         if( _calo->_roHalfTrans > 0.5*_calo->_crystalHW ) 
           {throw cet::exception("CaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}
        
      }
  }


}

