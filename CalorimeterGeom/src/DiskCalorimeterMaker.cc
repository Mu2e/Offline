//
// Make a Calorimeter.
//
// $Id: DiskCalorimeterMaker.cc,v 1.11 2014/02/25 01:09:42 echenard Exp $
// $Author: echenard $
// $Date: 2014/02/25 01:09:42 $

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
	_calo->_wrapperThickness      = config.getDouble("calorimeter.crystalWrapperThickness"); 
	_calo->_shellThickness        = config.getDouble("calorimeter.crystalShellThickness");

	_calo->_enveloppeInRadius     = config.getDouble("calorimeter.caloMotherInRadius"); 
	_calo->_enveloppeOutRadius    = config.getDouble("calorimeter.caloMotherOutRadius"); 
        _calo->_enveloppeZ0           = config.getDouble("calorimeter.caloMotherZ0"); 
        _calo->_enveloppeZ1           = config.getDouble("calorimeter.caloMotherZ1"); 

	_calo->_nROPerCrystal         = config.getInt("calorimeter.crystalReadoutChannelCount");
	_calo->_roHalfTrans           = config.getDouble("calorimeter.crystalReadoutHalfTrans");
	_calo->_roHalfThickness       = config.getDouble("calorimeter.crystalReadoutHalfThickness");
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

	_calo->_nPipes                = config.getInt("calorimeter.nPipes",0);    
	_calo->_pipeRadius            = config.getDouble("calorimeter.pipeRadius",5); 
	_calo->_pipeThickness         = config.getDouble("calorimeter.pipeThickness",0.5);   
	config.getVectorDouble("calorimeter.pipeTorRadius", _calo->_pipeTorRadius, _calo->_nPipes);
	
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

	// Create vanes
	MakeDisks();
 
 

  }


  DiskCalorimeterMaker::~DiskCalorimeterMaker() {}


  void DiskCalorimeterMaker::MakeDisks(void)
  {

      double diskHalfZLength    = _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness + _calo->_pipeRadius;
      double crystalCellRadius  = _calo->_crystalHalfTrans + _calo->_wrapperThickness + _calo->_shellThickness;      
      CLHEP::Hep3Vector crystalShiftInDisk(0,0,-_calo->_roHalfThickness);
      
      _calo->_nCrystalTot = 0;
            
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
	  _calo->_nCrystalTot += thisDisk->nCrystals();

  	  if (_verbosityLevel) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
        	                        <<" (X,Y,Z)="<<thisDisk->origin()<<"  local_(X,Y,Z)="<<thisDisk->originLocal()<<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;
          
	  if (_verbosityLevel > 1)
	  {
	      double espa          = thisDisk->estimateEmptySpace();
              double diskVolume    = 3.1415926*(thisDisk->outerRadius()*thisDisk->outerRadius()-thisDisk->innerRadius()*thisDisk->innerRadius());
              double crystalVolume = 3.4641016*crystalCellRadius*crystalCellRadius*thisDisk->nCrystals();

              std::cout<<"Estimated empty space between the disks and the crystals "<<std::endl;
	      std::cout<<"Inner edge and crystals = "<<espa<<std::endl;
              std::cout<<"Outer edge and crystals = "<<diskVolume-crystalVolume-espa<<" "<<std::endl;
	  }   
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
        
	if ( (_calo->_diskOuterRadius[i]+_calo->_caseThickness) > _calo->_enveloppeOutRadius) 
                    {throw cet::exception("DiskCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother for index="<<i<<".\n";} 

	if ( (_calo->_diskInnerRadius[i]-_calo->_caseThickness) < _calo->_enveloppeInRadius) 
                    {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius smaller than calorimeter mother for index="<<i<<".\n";} 

      }  
      
      double diskLength   = 2.0*( _calo->_crystalHalfLength + _calo->_roHalfThickness + _calo->_wrapperThickness + _calo->_caseThickness + _calo->_pipeRadius);
      double calozEnd     = _calo->_origin.z() + diskLength + _calo->_diskSeparation[_calo->_nSections-1]; 
      double calozBegin   = _calo->_origin.z();

      if (calozBegin < (_calo->_enveloppeZ0-0.1) || calozBegin > _calo->_enveloppeZ1) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.calorimeterZFront outside calorimeter mother (need 1mm margin for virtual detectors).\n";}  
      if (calozEnd   > (_calo->_enveloppeZ1-0.1))                       
          {throw cet::exception("DiskCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother (need 1mm margin for virtual detectors).\n";}  

  
      for (unsigned int i=0;i<_calo->_nPipes;++i) 
      {      
        if ( (_calo->_pipeTorRadius[i]- _calo->_pipeRadius) <  _calo->_enveloppeInRadius)           
	  {throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";}  

        if ( (_calo->_pipeTorRadius[i]+ _calo->_pipeRadius) >  _calo->_enveloppeOutRadius)           
	  {throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";}        
      }
      
 
  
  
  }


}//end mu2e namespace
