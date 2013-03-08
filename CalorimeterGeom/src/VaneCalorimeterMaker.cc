//
// Make a Vane Calorimeter.
//
// $Id: VaneCalorimeterMaker.cc,v 1.4 2013/03/08 01:22:31 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:31 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>


// Mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

// Framework include files
#include "cetlib/exception.h"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"


// other includes
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"





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

	_calo = std::auto_ptr<VaneCalorimeter>(new VaneCalorimeter());

	_calo->_nVane                = config.getInt   ("calorimeter.numberOfVanes");
	_calo->_crystalHW            = config.getDouble("calorimeter.crystalHalfTrans");
	_calo->_crystalHL            = config.getDouble("calorimeter.crystalHalfLong");
	_calo->_nCrystalR            = config.getInt   ("calorimeter.nCrystalRSlices");
	_calo->_nCrystalZ            = config.getInt   ("calorimeter.nCrystalZSlices");
	_calo->_rMin                 = config.getDouble("calorimeter.rInscribed");
	_calo->_rMax                 = _calo->_rMin+_calo->_nCrystalR*_calo->_crystalHW *2.0;



	_calo->_shieldHalfThickness  = config.getDouble("calorimeter.shieldHalfThickness");
	_calo->_neutronAbsorberHalfThickness  = config.getDouble("calorimeter.neutronAbsorberHalfThickness");

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
	_calo->_rMax                 = _calo->_rMin + 2.0*(_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR;

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

	// Local reference frame of the vane:
	// Z is along Mu2e z
	// Y is radial towards large radius
	// X is transverse to vane plane
    
        double crystalFullWidth = _calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness;
        CLHEP::Hep3Vector crystalShift(_calo->_roHalfThickness,0,0);

	double dX     =  _calo->_crystalHL + _calo->_wrapperThickness + _calo->_roHalfThickness;
	double dR     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR;
	double dZ     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalZ;
	double radius = _calo->_rMin + dR;
	double dphi   = 2*CLHEP::pi/_calo->_nVane;




	for( int i=0; i<_calo->_nVane; ++i ) {

          double phi = -CLHEP::pi + i*dphi;

          _calo->_vanes.push_back(Vane(i,_calo->_rMin,_calo->_nCrystalR,_calo->_nCrystalZ, crystalFullWidth, crystalShift));
          Vane& thisVane = _calo->_vanes.back();
          

          thisVane.setSize(        CLHEP::Hep3Vector(dX,dR,dZ) );
          thisVane.setOrigin(      CLHEP::Hep3Vector( _calo->_origin.x()+radius*cos(phi),_calo->_origin.y()+radius*sin(phi),_calo->_origin.z() ) );
          thisVane.setOriginLocal( CLHEP::Hep3Vector(radius*cos(phi),radius*sin(phi),0) );
          thisVane.setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(-CLHEP::pi/2 - i*dphi) );


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

