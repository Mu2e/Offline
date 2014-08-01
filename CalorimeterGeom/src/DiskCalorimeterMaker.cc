//
// Make a Calorimeter.
//
// $Id: DiskCalorimeterMaker.cc,v 1.13 2014/08/01 21:54:46 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 21:54:46 $

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
#include <memory>


// Mu2e includes
#include "cetlib/exception.h"
#include "CalorimeterGeom/inc/DiskCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
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

	_calo->_nSections  = config.getInt("calorimeter.numberOfDisks");      
        config.getVectorDouble("calorimeter.diskInnerRadius",  _calo->_diskInnerRadius, _calo->_nSections);
        config.getVectorDouble("calorimeter.diskOuterRadius",  _calo->_diskOuterRadius, _calo->_nSections);
        config.getVectorDouble("calorimeter.diskRotationAngle",_calo->_diskRotAngle,    _calo->_nSections);
        config.getVectorDouble("calorimeter.diskSeparation",   _calo->_diskSeparation,  _calo->_nSections);
	

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


	//THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
        double xOrigin                = -config.getDouble("mu2e.solenoidOffset");
        double zOrigin                = config.getDouble("calorimeter.calorimeterZFront",11750);
	_calo->_origin                = CLHEP::Hep3Vector(xOrigin,0,zOrigin);

        //standard formula to get the volume of the hexagon
	_calo->_caloGeomInfo.crystalVolume(6.9282032*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfLength());
		
        
	_verbosityLevel = config.getInt("calorimeter.verbosityLevel",0);

             
	//make sure above information is consistent
	CheckIt();

	// Create Disks
	MakeDisks();
 

  }


  DiskCalorimeterMaker::~DiskCalorimeterMaker() {}




  void DiskCalorimeterMaker::MakeDisks(void)
  {

      double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
      double crystalHalfTrans   = _calo->_caloGeomInfo.crystalHalfTrans();
      double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
      double caseThickness      = _calo->_caloGeomInfo.caseThickness();
      double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();
      double shellThickness     = _calo->_caloGeomInfo.shellThickness();
      double pipeRadius         = _calo->_caloGeomInfo.pipeRadius();
      
      
      double diskHalfZLength    = crystalHalfLength + roHalfThickness  + wrapperThickness + caseThickness + pipeRadius;
      double crystalCellRadius  = crystalHalfTrans  + wrapperThickness + shellThickness;      
      CLHEP::Hep3Vector crystalShiftInDisk(0,0,-roHalfThickness);
      


      for (unsigned int idisk=0; idisk<_calo->_nSections; ++idisk){			 
	 
	      CLHEP::Hep3Vector originLocal(0, 0, diskHalfZLength + _calo->_diskSeparation[idisk]);

	      double dR1    = _calo->_diskInnerRadius[idisk] - caseThickness;
	      double dR2    = _calo->_diskOuterRadius[idisk] + caseThickness;
	      double dZ     = 2.0*diskHalfZLength;

	      std::shared_ptr<Disk> thisDisk( new Disk(idisk,_calo->_diskInnerRadius[idisk], _calo->_diskOuterRadius[idisk], 2.0*crystalCellRadius, crystalShiftInDisk) );	 
	      _calo->_sections.push_back(thisDisk);

	      thisDisk->setSize(        CLHEP::Hep3Vector(dR1,dR2,dZ) );
              thisDisk->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(_calo->_diskRotAngle[idisk]) );
              thisDisk->setOriginLocal( originLocal );
              thisDisk->setOrigin(      originLocal + _calo->origin() );

  	      
	      //COORDINATE OF CENTER OF VOLUME CONTAINING THE CRYSTALS ONLY (NO READOUT,..) W.R.T CENTER OF OUTMOST DISK VOLUME
	      // outmost disk volume = volume ___originLocal vector____ points to, i.e outermost disk volume 
              CLHEP::Hep3Vector(0,0, pipeRadius - crystalHalfLength);
	      thisDisk->setCrystalShift( CLHEP::Hep3Vector(0,0, pipeRadius - crystalHalfLength) );
	      

	      //fill the full Crystal List / CaloSectionId (direct access to speed up computations) 
              int crystalOffset = _calo->_fullCrystalList.size();
	      for (int icry=0;icry<thisDisk->nCrystals();++icry)
	      {
	         
		 Crystal& thisCrystal = thisDisk->crystal(icry);
		 _calo->_fullCrystalList.push_back(&thisCrystal);
		 _calo->_crystalSectionId.push_back(idisk);
	      
	         //precompute the neighbors in the global frame
		 thisCrystal.setNeighborsLevel1(_calo->neighborsByLevel(icry+crystalOffset,1));
		 thisCrystal.setNeighborsLevel2(_calo->neighborsByLevel(icry+crystalOffset,2));
		 thisCrystal.setNeighborsLevel3(_calo->neighborsByLevel(icry+crystalOffset,3));
                 		 
                 //calculate the crystal position in the mu2e frame (aka global frame), taken from BaseCalorimeter.cc
		 CLHEP::Hep3Vector globalPosition = thisDisk->origin() + thisDisk->inverseRotation()*(thisCrystal.localPosition() + thisDisk->crystalShift()); 		 
		 thisCrystal.setPosition(globalPosition);

	      }	 


	      if (_verbosityLevel) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
        	                            <<" (X,Y,Z)="<<thisDisk->origin()<<"  local_(X,Y,Z)="<<thisDisk->originLocal()
					    <<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;          
	      
	      if (_verbosityLevel > 1)
	      {
		  double espace          = thisDisk->estimateEmptySpace();
        	  double diskVolume    = 3.1415926*(thisDisk->outerRadius()*thisDisk->outerRadius()-thisDisk->innerRadius()*thisDisk->innerRadius());
        	  double crystalVolume = 3.4641016*crystalCellRadius*crystalCellRadius*thisDisk->nCrystals();

        	  std::cout<<"Estimated empty space between the disks and the crystals "<<std::endl;
		  std::cout<<"Inner edge and crystals = "<<espace<<std::endl;
        	  std::cout<<"Outer edge and crystals = "<<diskVolume-crystalVolume-espace<<" "<<std::endl;
	      }   
            
      }
      
  }

 

  void DiskCalorimeterMaker::CheckIt(void)
  {
      
      int nROPerCrystal         = _calo->_caloGeomInfo.nROPerCrystal();
      double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
      double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
      double roHalfTrans        = _calo->_caloGeomInfo.roHalfTrans();
      double caseThickness      = _calo->_caloGeomInfo.caseThickness();
      double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();
      double pipeRadius         = _calo->_caloGeomInfo.pipeRadius();
 
      
      if( ! (nROPerCrystal ==1 || nROPerCrystal ==2 || nROPerCrystal ==4) ) 
          throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";      

      // Check size of readouts      
      if( nROPerCrystal==1 )
      {          
	  if (roHalfTrans > _calo->_caloGeomInfo.crystalHalfTrans() ) 
              throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";
        
      } else {
          if (roHalfTrans > 0.5*_calo->caloGeomInfo().crystalHalfTrans()) 
              throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";
      }
      
      
      
      //check enveloppe dimensions
      for (unsigned int i=0;i<_calo->_nSections;++i)
      {
        if (_calo->_diskInnerRadius[i] > _calo->_diskOuterRadius[i]) 
            throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius for index="<<i<<".\n";
        
	if ( (_calo->_diskOuterRadius[i] + caseThickness) > _calo->_caloGeomInfo.enveloppeOutRadius()) 
                    throw cet::exception("DiskCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother for index="<<i<<".\n";

	if ( (_calo->_diskInnerRadius[i] - caseThickness) < _calo->_caloGeomInfo.enveloppeInRadius()) 
                    {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius smaller than calorimeter mother for index="<<i<<".\n";} 

      }  
      
      
      //check disk length and envelope
      double diskLength   = 2.0*( crystalHalfLength + roHalfThickness + wrapperThickness + caseThickness + pipeRadius);
      double calozEnd     = _calo->_origin.z() + diskLength + _calo->_diskSeparation[_calo->_nSections-1]; 
      double calozBegin   = _calo->_origin.z();

      if (calozBegin < (_calo->_caloGeomInfo.enveloppeZ0()-0.1) || calozBegin > _calo->_caloGeomInfo.enveloppeZ1()) 
          throw cet::exception("DiskCaloGeom") << "calorimeter.calorimeterZFront outside calorimeter mother (need 1mm margin for virtual detectors).\n";
      
      if (calozEnd   > (_calo->_caloGeomInfo.enveloppeZ1()-0.1))                       
          throw cet::exception("DiskCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother (need 1mm margin for virtual detectors).\n";

  
      
      //look at pipes
      for (int i=0;i<_calo->caloGeomInfo().nPipes();++i) 
      {      
        if ( (_calo->_caloGeomInfo.pipeTorRadius().at(i)- _calo->_caloGeomInfo.pipeRadius()) <   _calo->_caloGeomInfo.enveloppeInRadius())           
	  throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";

        if ( (_calo->_caloGeomInfo.pipeTorRadius().at(i)+ _calo->_caloGeomInfo.pipeRadius()) >  _calo->_caloGeomInfo.enveloppeOutRadius())           
	  throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";        
      }
   
  }


}//end mu2e namespace
