//
// Make a Vane Calorimeter.
//
// $Id: VaneCalorimeterMaker.cc,v 1.12 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>
#include <memory>


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



// Vane geometry:
//
//  crystals + readout are places in the wrapping material
//  optional shell around the wrapping material, only on the sides of the crystal (no shell on front/back faces)
//  wrapper/shell are placed in the crystal box
//  shield / neutron absorber are in front of the crystal box
//  absorber + crystal box are placed inside the casing
//
//  cross section in z: casing - shield - absorber - wrapper - crystal - wrapper - crystal -... - wrapper - casing
//  cross section in height : casing - wrapper - crystal - readout - wrapper - casing
//  cross section in radius : casing - wrapper - crystal - wrapper - crystal - ...  - wrapper - casing

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

	_calo = std::unique_ptr<VaneCalorimeter>(new VaneCalorimeter());

	_calo->_nSections             = config.getInt("calorimeter.numberOfVanes");      
	_calo->_nCrystalR             = config.getInt   ("calorimeter.nCrystalRSlices");
	_calo->_nCrystalZ             = config.getInt   ("calorimeter.nCrystalZSlices");
	_calo->_shieldHalfThickness   = config.getDouble("calorimeter.shieldHalfThickness");
	_calo->_absorberHalfThickness = config.getDouble("calorimeter.neutronAbsorberHalfThickness");

	//Fill the Common Calo Data
	_calo->_caloGeomInfo.nROPerCrystal(      config.getInt("calorimeter.crystalReadoutChannelCount"));
	_calo->_caloGeomInfo.crystalHalfTrans(   config.getDouble("calorimeter.crystalHalfTrans") );
	_calo->_caloGeomInfo.crystalHalfLength(  config.getDouble("calorimeter.crystalHalfLong") );
	_calo->_caloGeomInfo.wrapperThickness(   config.getDouble("calorimeter.crystalWrapperThickness") );
	_calo->_caloGeomInfo.roHalfTrans(         config.getDouble("calorimeter.crystalReadoutHalfTrans") );
	_calo->_caloGeomInfo.roHalfThickness(    config.getDouble("calorimeter.crystalReadoutHalfThickness") );
	_calo->_caloGeomInfo.shellThickness(     config.getDouble("calorimeter.crystalShellThickness") );
	_calo->_caloGeomInfo.caseThickness(      config.getDouble("calorimeter.caseThickness") );
	_calo->_caloGeomInfo.enveloppeInRadius(  config.getDouble("calorimeter.caloMotherInRadius") );
	_calo->_caloGeomInfo.enveloppeOutRadius( config.getDouble("calorimeter.caloMotherOutRadius") );
	_calo->_caloGeomInfo.enveloppeZ0(        config.getDouble("calorimeter.caloMotherZ0") );
	_calo->_caloGeomInfo.enveloppeZ1(        config.getDouble("calorimeter.caloMotherZ1") );

	_calo->_caloGeomInfo.apdMeanNoise(       config.getDouble("calorimeter.meanNoiseAPD", 0.0) );
	_calo->_caloGeomInfo.apdSigmaNoise(      config.getDouble("calorimeter.sigmaNoiseAPD", 0.03) );
	_calo->_caloGeomInfo.lysoLightYield(     config.getDouble("calorimeter.lysoLightYield", 2000.0) );
	_calo->_caloGeomInfo.apdQuantumEff(      config.getDouble("calorimeter.quantumEffAPD", 0.68) );
	_calo->_caloGeomInfo.apdCollectEff(      config.getDouble("calorimeter.lightCollectEffAPD", 0.11));
	_calo->_caloGeomInfo.nonUniformity(      config.getDouble("calorimeter.crystalNonUniformity",0.0) );
	_calo->_caloGeomInfo.timeGap(            config.getDouble("calorimeter.timeGap",100.0) );
	_calo->_caloGeomInfo.electronEdep(       config.getDouble("calorimeter.electronDepositionAPD",1000.0) );
	_calo->_caloGeomInfo.electronEmin(       config.getDouble("calorimeter.electronMinEnergyAPD",0.1) );
		
	_calo->_caloGeomInfo.nPipes(             config.getInt("calorimeter.nPipes",0));
	_calo->_caloGeomInfo.pipeRadius(         config.getDouble("calorimeter.pipeRadius",5) );
	_calo->_caloGeomInfo.pipeThickness(      config.getDouble("calorimeter.pipeThickness",0.5) );
	
	std::vector<double> temp;
	config.getVectorDouble("calorimeter.pipeTorRadius", temp, _calo->_caloGeomInfo.nPipes());
	_calo->_caloGeomInfo.pipeTorRadius( temp );


        double crystalFullTrans    = _calo->_caloGeomInfo.crystalHalfTrans() + _calo->_caloGeomInfo.wrapperThickness() + _calo->_caloGeomInfo.shellThickness();
	_calo->_rMin               = config.getDouble("calorimeter.rInscribed");
	_calo->_rMax               = _calo->_rMin + 2.0*crystalFullTrans*_calo->_nCrystalR;


	//THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
        double xOrigin                = -config.getDouble("mu2e.solenoidOffset");
        double zOrigin                = config.getDouble("calorimeter.calorimeterZFront",11750);
	_calo->_origin                = CLHEP::Hep3Vector(xOrigin,0,zOrigin);

 
         //standard formula to get the volume of the rectangle
	_calo->_caloGeomInfo.crystalVolume(8*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfLength());

 	
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

          
	  double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
	  double crystalHalfTrans   = _calo->_caloGeomInfo.crystalHalfTrans();
	  double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
	  double caseThickness      = _calo->_caloGeomInfo.caseThickness();
	  double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();
	  double shellThickness     = _calo->_caloGeomInfo.shellThickness();
	  
	  double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
	  double crystalFullTrans   =  crystalHalfTrans + wrapperThickness + shellThickness;
	  double caloHalfLength     =  crystalFullTrans*_calo->_nCrystalZ;
          


	  double dX     = crystalHalfLength + wrapperThickness + roHalfThickness + caseThickness;
	  double dR     = crystalFullTrans * _calo->_nCrystalR + caseThickness;
	  double dZ     = crystalFullTrans * _calo->_nCrystalZ + absorberHalfLength + caseThickness;
	  double radius = _calo->_rMin + dR - caseThickness;
	  double dphi   = 2*CLHEP::pi/_calo->_nSections;
          CLHEP::Hep3Vector crystalShiftInDisk(-roHalfThickness,0,0);

	  for (unsigned int i=0; i<_calo->_nSections; ++i ) 
	  {
              double phi = -CLHEP::pi + i*dphi;

	      std::shared_ptr<Vane> thisVane( new Vane(i,_calo->_rMin,_calo->_nCrystalR,_calo->_nCrystalZ, crystalFullTrans, crystalShiftInDisk) );	 
	      _calo->_sections.push_back(thisVane);

	      CLHEP::Hep3Vector localOrigin(radius*cos(phi),radius*sin(phi),absorberHalfLength + caloHalfLength + caseThickness);

              thisVane->setSize(        CLHEP::Hep3Vector(dX,dR,dZ) );
              thisVane->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(CLHEP::pi/2 - i*dphi) );
              thisVane->setOriginLocal( localOrigin );
              thisVane->setOrigin(      localOrigin + _calo->origin() );             
	      thisVane->setCrystalShift( CLHEP::Hep3Vector(-_calo->_caloGeomInfo.crystalHalfLength(), 0 ,absorberHalfLength) );

	      //fill the full Crystal List (direct access to crystal from calorimeter as requested from users) 
              int crystalOffset = _calo->_fullCrystalList.size();
              for (int icry=0;icry<thisVane->nCrystals();++icry)
	      {
		 Crystal& thisCrystal = thisVane->crystal(icry);
		 _calo->_fullCrystalList.push_back(&thisCrystal);
		 _calo->_crystalSectionId.push_back(i);
	      
	         //precompute the neighbors in the global frame
		 thisCrystal.setNeighborsLevel1(_calo->neighborsByLevel(icry+crystalOffset,1));
		 thisCrystal.setNeighborsLevel2(_calo->neighborsByLevel(icry+crystalOffset,2));
		 thisCrystal.setNeighborsLevel3(_calo->neighborsByLevel(icry+crystalOffset,3));
                 thisCrystal.setPosition(_calo->crystalOrigin(icry));

                 //calculate the crystal position in the mu2e frame (aka global frame), taken from BaseCalorimeter.cc
		 CLHEP::Hep3Vector globalPosition = thisVane->origin() + thisVane->inverseRotation()*(thisCrystal.localPosition() + thisVane->crystalShift()); 		 
		 thisCrystal.setPosition(globalPosition);
	      }	 
	  }

      }



      void VaneCalorimeterMaker::CheckIt(void)
      {

  	   
	    //check that calorimeter fits in the mother volume
  	    double crystalFullTrans   =  _calo->_caloGeomInfo.crystalHalfTrans() + _calo->_caloGeomInfo.wrapperThickness() + _calo->_caloGeomInfo.shellThickness();
 	    double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
	    double dR                 = crystalFullTrans * _calo->_nCrystalR + _calo->_caloGeomInfo.caseThickness();
	    double dZ                 = crystalFullTrans * _calo->_nCrystalZ + absorberHalfLength + _calo->_caloGeomInfo.caseThickness();
            double calozBegin         = _calo->_origin.z();
            double calozEnd           = _calo->_origin.z() + 2*dZ;


	    if ( (_calo->_rMin + 2*dR) > _calo->_caloGeomInfo.enveloppeOutRadius()) 
                     {throw cet::exception("VaneCaloGeom") << "calorimeter outer radius larger than calorimeter mother \n";} 

	    if (  _calo->_rMin < _calo->_caloGeomInfo.enveloppeInRadius()) 
                     {throw cet::exception("VaneCaloGeom") << "calorimeter inner radius smaller than calorimeter mother \n";} 

	    if (calozBegin <  _calo->_caloGeomInfo.enveloppeZ0() || calozBegin >  _calo->_caloGeomInfo.enveloppeZ1()) 
                     {throw cet::exception("VaneCaloGeom") << "calorimeter.calorimeterZFront   outside calorimeter mother.\n";}  

	    if (calozEnd   > _calo->_caloGeomInfo.enveloppeZ1())                       
                     {throw cet::exception("VaneCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother.\n";}  


	    // Check number of readouts
	    int nRO    = _calo->caloGeomInfo().nROPerCrystal();
 	    if( ! (nRO==1 || nRO==2 || nRO==4) ) 
              {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}


	    // Check size of readouts
	    if( nRO==1 ) {
               if( _calo->_caloGeomInfo.roHalfTrans() > _calo->_caloGeomInfo.crystalHalfTrans() ) 
        	 {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHalfTrans.\n";}

	    } else {
               if( _calo->_caloGeomInfo.roHalfTrans()> 0.5*_calo->_caloGeomInfo.crystalHalfTrans() ) 
        	 {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}

	    }
            
      }


}

