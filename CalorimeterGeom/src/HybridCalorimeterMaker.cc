//
// Make a Calorimeter.
//
// $Id: HybridCalorimeterMaker.cc,v 1.2 2014/02/25 01:09:42 echenard Exp $
// $Author: echenard $
// $Date: 2014/02/25 01:09:42 $

// original authors Julie Managan and Robert Bernstein
// quite a few changes by G. Pezzullo

// Disk geometry
//
//  crystals + readout are places in the wrapping material
//  optional shell around the wrapping material, only on the sides of the crystal (no shell on front/back faces)
//  wrapper/shell are placed in the disk
//  disk is placed into the casing (i.e. a bigger disk)
//  outermost layer is the casing
//
//  cross section in z: casing - wrapper - crystal - readout - wrapper - casing
//  cross section in radius: casing - wrapper - crystal - wrapper - crystal - ... - wrapper - casing
 
//Barrel Geometry
//
//
//
//

// C++ includes
#include <math.h>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "cetlib/exception.h"
#include "CalorimeterGeom/inc/HybridCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Barrel.hh"
#include "CalorimeterGeom/inc/HybridCalorimeter.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"

// other includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"




namespace mu2e{


  HybridCalorimeterMaker::HybridCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
  {

    _calo = std::unique_ptr<HybridCalorimeter>(new HybridCalorimeter());

    _calo->_nSections             = 2.0;//config.getInt("calorimeter.numberOfDisks");      
    _calo->_caseThickness         = config.getDouble("calorimeter.caseThickness");
    _calo->_diskInnerRadius       = config.getDouble("calorimeter.diskInnerRadius");
    _calo->_diskOuterRadius       = config.getDouble("calorimeter.diskOuterRadius");
    _calo->_diskRotAngle          = config.getDouble("calorimeter.diskRotationAngle");

    //wheel parameters
    _calo->_nWheels               = config.getDouble("calorimeter.nWheels");
    _calo->_nCrystalWheel         = config.getDouble("calorimeter.nCrystalWheel");
    _calo->_barrelInnerRadius     = config.getDouble("calorimeter.barrelInnerRadius");
    _calo->_diskToBarrelSeparation= config.getDouble("calorimeter.diskToBarrelSeparation");
    _calo->_barrelRotAngle        = config.getDouble("calorimeter.barrelRotAngle");
    
    _calo->_hexCrystalHalfTrans   = config.getDouble("calorimeter.hexCrystalHalfTrans");
    _calo->_hexCrystalHalfLength  = config.getDouble("calorimeter.hexCrystalHalfLong");
    _calo->_hexCrystalVolume      = 3.4641016*_calo->_hexCrystalHalfTrans*_calo->_hexCrystalHalfTrans*2.0*_calo->_hexCrystalHalfLength;//FIX ME
    _calo->_crystalHalfTrans      = config.getDouble("calorimeter.crystalHalfTrans");
    _calo->_crystalHalfLength     = config.getDouble("calorimeter.crystalHalfLong");
    double tanTheta = std::tan(CLHEP::pi / _calo->_nCrystalWheel);
    double Bmin  = _calo->_barrelInnerRadius*tanTheta*2.0;
    double Bmax  = (_calo->_barrelInnerRadius + 2.0*_calo->_crystalHalfLength)*tanTheta*2.0;
    _calo->_barrelCrystalVolume   = (Bmin + Bmax)*_calo->_crystalHalfLength*_calo->_nWheels*_calo->_nCrystalWheel;
    _calo->_wrapperThickness      = config.getDouble("calorimeter.crystalWrapperThickness",0.0); 
    _calo->_shellThickness        = config.getDouble("calorimeter.crystalShellThickness",0.0);
    _calo->_enveloppeInRadius     = config.getDouble("calorimeter.caloMotherInRadius",0); 
    _calo->_enveloppeOutRadius    = config.getDouble("calorimeter.caloMotherOutRadius",765); 
    _calo->_enveloppeZ0           = config.getDouble("calorimeter.caloMotherZ0",11740); 
    _calo->_enveloppeZ1           = config.getDouble("calorimeter.caloMotherZ1",13910); 

    _calo->_nROPerCrystal         = config.getInt("calorimeter.crystalReadoutChannelCount");
    _calo->_roHalfTrans           = config.getDouble("calorimeter.crystalReadoutHalfTrans");
    _calo->_roHalfThickness       = config.getDouble("calorimeter.crystalReadoutHalfThickness");

    _calo->_barrelOuterRadius     = _calo->_barrelInnerRadius + 2.0*( _calo->_crystalHalfLength + _calo->_wrapperThickness + _calo->_shellThickness + _calo->_caseThickness + _calo->_roHalfThickness);

    _calo->_pipeRadius            = config.getDouble("calorimeter.pipeRadius",5); 
    _calo->_pipeThickness         = config.getDouble("calorimeter.pipeThickness",0.25); 

    _calo->_nonUniformity         = config.getDouble("calorimeter.crystalNonUniformity",0.0);
    _calo->_timeGap               = config.getDouble("calorimeter.timeGap",100.0);
    _calo->_electronEdep          = config.getDouble("calorimeter.electronDepositionAPD",1000.0);
    _calo->_electronEmin          = config.getDouble("calorimeter.electronMinEnergyAPD",0.1);

    _calo->_apdMeanNoise          = config.getDouble("calorimeter.meanNoiseAPD", 0.0);
    _calo->_apdSigmaNoise         = config.getDouble("calorimeter.sigmaNoiseAPD", 0.03);
    _calo->_lysoLightYield        = config.getDouble("calorimeter.lysoLightYield", 2000.0);
    _calo->_apdQuantumEff         = config.getDouble("calorimeter.quantumEffAPD", 0.68);
    _calo->_lightCollectEffAPD    = config.getDouble("calorimeter.lightCollectEffAPD", 0.11);

    _verbosityLevel               = config.getInt("calorimeter.verbosityLevel",0);

     

    //THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
    double xOrigin                = -config.getDouble("mu2e.solenoidOffset");
    double zOrigin                = config.getDouble("calorimeter.calorimeterZFront",11750);
    _calo->_origin                = CLHEP::Hep3Vector(xOrigin,0,zOrigin);
     
     
    //COORDINATE OF CENTER OF VOLUME CONTAINING THE CRYSTALS ONLY (NO READOUT,..) W.R.T CENTER OF OUTMOST DISK VOLUME
    // outmost disk volume = volume ___originLocal vector____ points to, i.e outermost disk volume 
    _calo->_crystalShift = CLHEP::Hep3Vector(0,0, _calo->_pipeRadius -_calo->_crystalHalfLength);
   
          
    //make sure above information is consistent
    CheckIt();

    // Create disk and barrel
    _calo->_nCrystalTot = 0;
    MakeDisk();
    MakeBarrel();

  }


  HybridCalorimeterMaker::~HybridCalorimeterMaker() {}


  void HybridCalorimeterMaker::MakeDisk(void)
  {

    double diskHalfZLength    = _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness + _calo->_pipeRadius;
    double crystalCellRadius  = _calo->_crystalHalfTrans + _calo->_wrapperThickness + _calo->_shellThickness;      
    CLHEP::Hep3Vector crystalShiftInDisk(0,0,-_calo->_roHalfThickness);
            
    CLHEP::Hep3Vector originLocal(0, 0, diskHalfZLength);

    double dR1    = _calo->_diskInnerRadius - _calo->_caseThickness;
    double dR2    = _calo->_diskOuterRadius + _calo->_caseThickness;
    double dZ     = 2.0*diskHalfZLength;

    boost::shared_ptr<Disk> thisDisk( new Disk(0,_calo->_diskInnerRadius, _calo->_diskOuterRadius, 2.0*crystalCellRadius, crystalShiftInDisk) );	 
    _calo->_sections.push_back(thisDisk);

    thisDisk->setSize(        CLHEP::Hep3Vector(dR1,dR2,dZ) );
    thisDisk->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(_calo->_diskRotAngle) );
    thisDisk->setOriginLocal( originLocal );
    thisDisk->setOrigin(      originLocal + _calo->origin() );
    _calo->_nCrystalTot += thisDisk->nCrystals();

    if (_verbosityLevel) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
				  <<" (X,Y,Z)="<<thisDisk->origin()<<"  local_(X,Y,Z)="<<thisDisk->originLocal()<<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;

      

  }

void HybridCalorimeterMaker::MakeBarrel(void)
  {

    
    double crystalCellRadius  = _calo->_crystalHalfTrans + _calo->_wrapperThickness + _calo->_shellThickness;     
    double barrelHalfZLength    = crystalCellRadius*_calo->_nWheels + 2.0*_calo->_caseThickness; 
    CLHEP::Hep3Vector crystalShiftInDisk(0,0,-_calo->_roHalfThickness);
            
    CLHEP::Hep3Vector originLocal(0, 0, barrelHalfZLength + _calo->_diskToBarrelSeparation);

    double dR1    = _calo->_barrelInnerRadius - _calo->_caseThickness;
    double dR2    = _calo->_barrelOuterRadius + _calo->_caseThickness;
    double dZ     = 2.0*barrelHalfZLength;

    boost::shared_ptr<Barrel> thisBarrel( new Barrel(1,_calo->_barrelInnerRadius, _calo->_barrelOuterRadius, 2.0*crystalCellRadius, _calo->_nWheels, _calo->_nCrystalWheel,  crystalShiftInDisk) );	 
    _calo->_sections.push_back(thisBarrel);

    thisBarrel->setSize(        CLHEP::Hep3Vector(dR1,dR2,dZ) );
    thisBarrel->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(_calo->_barrelRotAngle) );
    thisBarrel->setOriginLocal( originLocal );
    thisBarrel->setOrigin(      originLocal + _calo->origin() );
    _calo->_nCrystalTot += thisBarrel->nCrystals();

    if (_verbosityLevel) std::cout<<"Constructed Barrel "<<thisBarrel->id()<<":  Rin="<<thisBarrel->innerRadius()<<"  Rout="<<thisBarrel->outerRadius()
				  <<" (X,Y,Z)="<<thisBarrel->origin()<<"  local_(X,Y,Z)="<<thisBarrel->originLocal()<<"  with "<<thisBarrel->nCrystals()<<" crystals"<<std::endl;

      

  }
 
 
  void HybridCalorimeterMaker::CheckIt(void)
  {
      
    if( ! (_calo->_nROPerCrystal ==1 || _calo->_nROPerCrystal ==2 || _calo->_nROPerCrystal ==4) ) 
      {throw cet::exception("HybridCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}      

    // Check size of readouts      
    if( _calo->_nROPerCrystal==1 ) {
      if(  _calo->_roHalfTrans > _calo->_crystalHalfTrans ) 
	{throw cet::exception("BarrelCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";}
        
    } else {
      if( _calo->_roHalfTrans > 0.5*_calo->_crystalHalfTrans) 
	{throw cet::exception("BarrelCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}
    }
      
    
    if (_calo->_diskInnerRadius > _calo->_diskOuterRadius) 
      {throw cet::exception("BarrelCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius\n";} 
    
    if ( (_calo->_diskOuterRadius+_calo->_caseThickness) > _calo->_enveloppeOutRadius) 
      {throw cet::exception("BarrelCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother\n";} 
    
    if ( (_calo->_diskInnerRadius-_calo->_caseThickness) < _calo->_enveloppeInRadius) 
      {throw cet::exception("BarrelCaloGeom") << "calorimeter.diskInnerRadius smaller than calorimeter mother\n";} 
    
    
    
    //double diskLength   = 2.0*( _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness + _calo->_pipeRadius);

    double barrelLength = 2.0*( _calo->_wrapperThickness + _calo->_caseThickness + _calo->_crystalHalfTrans)*_calo->_nWheels;
    double calozEnd     = _calo->_origin.z() + barrelLength + _calo->_diskToBarrelSeparation;
    double calozBegin   = _calo->_origin.z();

    if (calozBegin < (_calo->_enveloppeZ0-0.1) || calozBegin > _calo->_enveloppeZ1) 
      {throw cet::exception("HybridCaloGeom") << "calorimeter.calorimeterZFront outside calorimeter mother (need 1mm margin for virtual detectors).\n";}  
    if (calozEnd   > (_calo->_enveloppeZ1-0.1))                       
      {throw cet::exception("HybridCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother (need 1mm margin for virtual detectors).\n";}  

  }


}//end mu2e namespace
