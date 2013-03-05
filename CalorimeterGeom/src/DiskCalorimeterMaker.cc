//
// Make a Calorimeter.
//
// $Id: DiskCalorimeterMaker.cc,v 1.3 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>

//
// Mu2e includes
#include "CalorimeterGeom/inc/DiskCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

// Framework include files
#include "cetlib/exception.h"
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


  DiskCalorimeterMaker::DiskCalorimeterMaker( SimpleConfig const& config, double solenoidOffset)
  {

        _calo = std::auto_ptr<DiskCalorimeter>(new DiskCalorimeter());

	_calo->_nDisks                = config.getInt("calorimeter.numberOfDisks");      
	_calo->_diskThickness         = config.getDouble("calorimeter.diskCaseThickness");
        config.getVectorDouble("calorimeter.diskInnerRadius",_calo->_diskInnerRadius,_calo->_nDisks);
        config.getVectorDouble("calorimeter.diskOuterRadius",_calo->_diskOuterRadius,_calo->_nDisks);
        config.getVectorDouble("calorimeter.diskRotationAngle",_calo->_diskRotAngle,_calo->_nDisks);
        config.getVectorDouble("calorimeter.diskSeparation",_calo->_diskSeparation,_calo->_nDisks);

	_calo->_crystalHalfTrans      = config.getDouble("calorimeter.crystalHalfTrans");
	_calo->_crystalDepth          = 2.0*config.getDouble("calorimeter.crystalHalfLong");  //half length in config file...       
	_calo->_wrapperThickness      = config.getDouble("calorimeter.crystalWrapperThickness",0.0); 
	_calo->_shellThickness        = config.getDouble("calorimeter.crystalShellThickness",0.0);


	_calo->_nROPerCrystal         = config.getInt("calorimeter.crystalReadoutChannelCount");
	_calo->_readOutHalfSize       = config.getDouble("calorimeter.crystalReadoutHalfTrans");
	_calo->_readOutHalfThickness  = config.getDouble("calorimeter.crystalReadoutHalfThickness");

	_calo->_nonUniformity         = config.getDouble("calorimeter.crystalNonUniformity",0.0);
	_calo->_timeGap               = config.getDouble("calorimeter.timeGap",100.0);
	_calo->_electronEdep          = config.getDouble("calorimeter.electronDepositionAPD",1000.0);
	_calo->_electronEmin          = config.getDouble("calorimeter.electronMinEnergyAPD",0.1);

        verbosityLevel                = config.getInt("calorimeter.verbosityLevel",0);


        //calorimeter center refers to center of first disk
        CLHEP::Hep3Vector center = config.getHep3Vector("calorimeter.calorimeterCenter");
        //in construct calorimeter, the calo is placed at center.x(),center.y(),center.z()-zOffset
	//zoffset is the position of toyDS3 in Mu2e coordinates. //this is set in geom_01.txt
	
        _calo->_origin = CLHEP::Hep3Vector(center.x()-solenoidOffset,center.y(),center.z());

	//make sure above information is consistent
	CheckIt();

	// Create vanes
	MakeDisks();


  }


  DiskCalorimeterMaker::~DiskCalorimeterMaker() {}


 
  void DiskCalorimeterMaker::MakeDisks(void)
  {

      double crystalFullWidth = _calo->_crystalHalfTrans + _calo->_wrapperThickness + _calo->_shellThickness;
      CLHEP::Hep3Vector crystalShift(0,0,-_calo->_readOutHalfThickness);

      for (unsigned int idisk=0; idisk<_calo->_nDisks; ++idisk)
      {			 
	 
	 double dR1    = _calo->_diskInnerRadius[idisk] - _calo->_diskThickness;
	 double dR2    = _calo->_diskOuterRadius[idisk] + _calo->_diskThickness;
	 double dZ     = _calo->_crystalDepth + 2.0*_calo->_readOutHalfThickness + 2.0*_calo->_diskThickness;
	 
         


	 _calo->_disks.push_back(Disk(idisk,_calo->_diskInnerRadius[idisk],_calo->_diskOuterRadius[idisk],_calo->_diskThickness,2.0*crystalFullWidth, crystalShift));
	 Disk &thisDisk = _calo->_disks.back();
	 	 
         thisDisk.setSize(        CLHEP::Hep3Vector(dR1,dR2,dZ) );
         thisDisk.setOrigin(      CLHEP::Hep3Vector(_calo->_origin.x(), _calo->_origin.y(), _calo->_origin.z()+_calo->_diskSeparation[idisk]) );
         thisDisk.setOriginLocal( CLHEP::Hep3Vector(0,0,_calo->_diskSeparation[idisk]) );
         thisDisk.setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(_calo->_diskRotAngle[idisk]) );
	 	 
  	 if ( verbosityLevel > 0) std::cout<<"Constructed Disk "<<thisDisk.id()<<":  Rin="<<thisDisk.innerRadius()<<"  Rout="<<thisDisk.outerRadius()
        	  <<" X="<<thisDisk.origin().x()
        	  <<" Y="<<thisDisk.origin().y()
        	  <<" Z="<<thisDisk.origin().z()
		  <<"  with "<<thisDisk.nCrystals()<<" crystals"<<std::endl;
         }

  }

 
 
  void DiskCalorimeterMaker::CheckIt(void)
  {


      if( ! (_calo->_nROPerCrystal ==1 || _calo->_nROPerCrystal ==2 || _calo->_nROPerCrystal ==4) ) 
        {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}      

      // Check size of readouts      
      if( _calo->_nROPerCrystal==1 ) {
        if(  _calo->_readOutHalfSize > _calo->_crystalHalfTrans ) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";}
        
      } else {
        if( _calo->_readOutHalfSize > 0.5*_calo->_crystalHalfTrans) 
          {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}
      }
      
      for (unsigned int i=0;i<_calo->_nDisks;++i) {
        if (_calo->_diskInnerRadius[i] > _calo->_diskOuterRadius[i]) 
            {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius for index="<<i<<".\n";} 
      }  


  }


}//end mu2e namespace
