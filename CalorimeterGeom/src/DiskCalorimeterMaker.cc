//
// Make a Calorimeter.
//
// $Id: DiskCalorimeterMaker.cc,v 1.6 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $

// original authors Julie Managan and Robert Bernstein
// quite a few changes by Bertrand Echenarrd

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
 

// C++ includes
#include <math.h>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "cetlib/exception.h"
#include "CalorimeterGeom/inc/DiskCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
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


  DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
  {

        _calo = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());

	_calo->_nSections             = config.getInt("calorimeter.numberOfDisks");      
	_calo->_caseThickness         = config.getDouble("calorimeter.caseThickness");
        config.getVectorDouble("calorimeter.diskInnerRadius",  _calo->_diskInnerRadius, _calo->_nSections);
        config.getVectorDouble("calorimeter.diskOuterRadius",  _calo->_diskOuterRadius, _calo->_nSections);
        config.getVectorDouble("calorimeter.diskRotationAngle",_calo->_diskRotAngle,    _calo->_nSections);
        config.getVectorDouble("calorimeter.diskSeparation",   _calo->_diskSeparation,  _calo->_nSections);

	_calo->_crystalHalfTrans      = config.getDouble("calorimeter.crystalHalfTrans");
	_calo->_crystalHalfLength     = config.getDouble("calorimeter.crystalHalfLong");
	_calo->_wrapperThickness      = config.getDouble("calorimeter.crystalWrapperThickness",0.0); 
	_calo->_shellThickness        = config.getDouble("calorimeter.crystalShellThickness",0.0);

	_calo->_enveloppeRadius       = config.getDouble("calorimeter.caloMotherRadius",850); 
        _calo->_enveloppeZ0           = config.getDouble("calorimeter.caloMotherZ0",11740); 
        _calo->_enveloppeZ1           = config.getDouble("calorimeter.caloMotherZ1",13910); 

	_calo->_nROPerCrystal         = config.getInt("calorimeter.crystalReadoutChannelCount");
	_calo->_roHalfTrans           = config.getDouble("calorimeter.crystalReadoutHalfTrans");
	_calo->_roHalfThickness       = config.getDouble("calorimeter.crystalReadoutHalfThickness");

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
        double zOrigin                = config.getDouble("calorimeter.calorimeterZOrigin",11740);
	_calo->_origin                = CLHEP::Hep3Vector(xOrigin,0,zOrigin);
     
     
	//COORDINATE OF VOLUME CONTAINING THE CRYSTALS ONLY (NO READOUT,..) W.R.T OUTSIDE DISK VOLUME
	// outside disk volume = volume ___originLocal vector____ points to, i.e outermost disk volume 
        _calo->_crystalShift = CLHEP::Hep3Vector(0,0,-_calo->_crystalHalfLength);
       
	//make sure above information is consistent
	CheckIt();

	// Create vanes
	MakeDisks();


  }


  DiskCalorimeterMaker::~DiskCalorimeterMaker() {}


  void DiskCalorimeterMaker::MakeDisks(void)
  {

      double diskHalfZLength    = _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness;
      double crystalCellRadius  = _calo->_crystalHalfTrans + _calo->_wrapperThickness + _calo->_shellThickness;      
      CLHEP::Hep3Vector crystalShiftInDisk(0,0,-_calo->_roHalfThickness);
      
            
      for (unsigned int idisk=0; idisk<_calo->_nSections; ++idisk)
      {			 
	 
	  CLHEP::Hep3Vector originLocal(0, 0, diskHalfZLength + _calo->_diskSeparation[idisk]);

	  double dR1    = _calo->_diskInnerRadius[idisk] - _calo->_caseThickness;
	  double dR2    = _calo->_diskOuterRadius[idisk] + _calo->_caseThickness;
	  double dZ     = 2.0*diskHalfZLength;

	  boost::shared_ptr<Disk> thisDisk( new Disk(idisk,_calo->_diskInnerRadius[idisk], _calo->_diskOuterRadius[idisk], 2.0*crystalCellRadius, crystalShiftInDisk) );	 
	  _calo->_sections.push_back(thisDisk);

	  thisDisk->setSize(        CLHEP::Hep3Vector(dR1,dR2,dZ) );
          thisDisk->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(_calo->_diskRotAngle[idisk]) );
          thisDisk->setOriginLocal( originLocal );
          thisDisk->setOrigin(      originLocal + _calo->origin() );

  	  if (_verbosityLevel) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
        	                        <<" (X,Y,Z)="<<thisDisk->origin()<<"  local_(X,Y,Z)="<<thisDisk->originLocal()<<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;

      }

  }

 
 
  void DiskCalorimeterMaker::CheckIt(void)
  {
      
      if( ! (_calo->_nROPerCrystal ==1 || _calo->_nROPerCrystal ==2 || _calo->_nROPerCrystal ==4) ) 
        {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}      

      // Check size of readouts      
      if( _calo->_nROPerCrystal==1 ) {
        if(  _calo->_roHalfTrans > _calo->_crystalHalfTrans ) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";}
        
      } else {
        if( _calo->_roHalfTrans > 0.5*_calo->_crystalHalfTrans) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}
      }
      
      for (unsigned int i=0;i<_calo->_nSections;++i) {
        if (_calo->_diskInnerRadius[i] > _calo->_diskOuterRadius[i]) 
            {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius for index="<<i<<".\n";} 
        
	if ( (_calo->_diskOuterRadius[i]+_calo->_caseThickness) > _calo->_enveloppeRadius) 
                    {throw cet::exception("DiskCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother for index="<<i<<".\n";} 

      }  
      
      double diskLength   = 2.0*( _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness);
      double calozEnd     = _calo->_origin.z() + diskLength + _calo->_diskSeparation[_calo->_nSections-1]; 
      double calozBegin   = _calo->_origin.z();

      if (calozBegin < _calo->_enveloppeZ0 || calozBegin > _calo->_enveloppeZ1) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.calorimeterZOrigin   outside calorimeter mother.\n";}  
      if (calozEnd   > _calo->_enveloppeZ1)                       
          {throw cet::exception("DiskCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother.\n";}  

  }


}//end mu2e namespace
