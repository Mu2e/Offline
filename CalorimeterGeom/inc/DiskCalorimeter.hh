#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// $Id: DiskCalorimeter.hh,v 1.3 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Hold all geometry and identifier information about
// a Disk Calorimeter. In order to insulate this class from
// knowledge of databases etc, this class can not know
// how to make itself.
//

// Look at Mu2eG4/inc/constructDiskCalorimeter.cc 
// for definition of geometry

// Original author B. Echenard


//C++ includes
#include <vector>

// Mu2e includes
#include "Calorimeter.hh"
#include "Disk.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {


    class DiskCalorimeter: public Calorimeter{

       friend class DiskCalorimeterMaker;


       public:


         DiskCalorimeter()  {}
         ~DiskCalorimeter() {}


	 unsigned int nDisks(void) const        {return _nDisks;}
	 Disk const& disk(int i) const       {return _disks.at(i);}
	 double diskSeparation(int i) const  {return _diskSeparation.at(i);}

	 double envelopeRmin(void) const ;
	 double envelopeRmax(void) const ;
	 double envelopeHalfLength(void) const ;

         unsigned int nROPerCrystal(void) const {return _nROPerCrystal;}
         double crystalVolume(void) const       {return 3.4641016*_crystalHalfTrans*_crystalHalfTrans*_crystalDepth;}
         unsigned int nRO(void) const;
         unsigned int nCrystal(void) const;


         int caloSectionId(int crystalId) const;
         int localCrystalId(int crystalId) const;

         //conversion of crystal <-> readout id, 
	 //readout_id = crystal_id*nRoPerCrystal ... crystal_id*nRoPerCrystal + nRoPerCrystal-1		 
	 int crystalByRO(int roid) const     {return (roid/_nROPerCrystal);}
	 int ROBaseByCrystal(int crystalId) const   {return (crystalId*_nROPerCrystal);}
	 
         std::vector<int> neighbors(int crystalId, int level=1) const; 


	 CLHEP::Hep3Vector const& origin(void) const  {return _origin;}

	 CLHEP::Hep3Vector crystalOrigin(int crystalId) const;
	 CLHEP::Hep3Vector localCrystalOrigin(int crystalId) const;
         CLHEP::Hep3Vector crystalAxis(int crystalId) const ;
	 CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;


	 double roHalfThickness(void) const     {return _readOutHalfThickness;}
	 double roHalfSize(void) const          {return _readOutHalfSize;}

	 double crysHalfTrans(void) const       {return _crystalHalfTrans;}
	 double crystalHalfLength(void) const           {return _crystalDepth;}
     double wrapperThickness(void) const    {return _wrapperThickness; }
	 double shellThickness(void) const      {return _shellThickness;}

	 double getNonuniformity(void) const    {return _nonUniformity; }
	 double getTimeGap(void) const          {return _timeGap; }
	 double getElectronEdep(void) const     {return _electronEdep; }
	 double getElectronEmin(void) const     {return _electronEmin; }




         	 
 	 	 



       protected:

	  unsigned int _nDisks;
	  double       _diskThickness;          
	  std::vector<double> _diskInnerRadius;
	  std::vector<double> _diskOuterRadius;
	  std::vector<double> _diskSeparation;
	  std::vector<double> _diskRotAngle;  


	  unsigned int _nROPerCrystal;
	  double _readOutHalfThickness;
	  double _readOutHalfSize;

	  double _crystalHalfTrans;
	  double _crystalDepth;
	  double _wrapperThickness;
	  double _shellThickness;

	  double _nonUniformity;
	  double _timeGap;
	  double _electronEdep; // energy deposition of charged particle crossing APD
	  double _electronEmin; // minimum energy deposition of charged particle crossing APD

	  std::vector<Disk>  _disks;
	  CLHEP::Hep3Vector  _origin;
	  	

     };

}    

#endif /* CalorimeterGeom_DiskCalorimeter_hh */
