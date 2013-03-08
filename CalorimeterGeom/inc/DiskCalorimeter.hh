#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// $Id: DiskCalorimeter.hh,v 1.4 2013/03/08 01:22:31 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:31 $
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



	 unsigned int nDisk() const          {return _nDisk;}
	 Disk const& disk(int i) const       {return _disks.at(i);}
	 double diskSeparation(int i) const  {return _diskSeparation.at(i);}

         unsigned int nROPerCrystal() const  {return _nROPerCrystal;}
         double crystalVolume() const        {return 3.4641016*_crystalHalfTrans*_crystalHalfTrans*_crystalDepth;}
         unsigned int nRO() const;
         unsigned int nCrystal() const;


         int caloSectionId(int crystalId) const;
         int localCrystalId(int crystalId) const;

         //conversion of crystal <-> readout id, 
	 //readout_id = crystal_id*nRoPerCrystal ... crystal_id*nRoPerCrystal + nRoPerCrystal-1		 
	 int crystalByRO(int roid) const            {return (roid/_nROPerCrystal);}
	 int ROBaseByCrystal(int crystalId) const   {return (crystalId*_nROPerCrystal);}
	 
         std::vector<int> neighbors(int crystalId, int level=1) const; 


	 CLHEP::Hep3Vector const& origin() const   {return _origin;}

	 CLHEP::Hep3Vector crystalOrigin(int crystalId) const;
	 CLHEP::Hep3Vector localCrystalOrigin(int crystalId) const;
         CLHEP::Hep3Vector crystalAxis(int crystalId) const ;
	 CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;

	 bool isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const ;       	 
	 bool isInsideDisk(int idisk, CLHEP::Hep3Vector const& pos) const ;       	 
	 int  crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const ;



	 double envelopeRmin() const ;
	 double envelopeRmax() const ;
	 double envelopeHalfLength() const ;

	 double roHalfThickness() const     {return _readOutHalfThickness;}
	 double roHalfSize() const          {return _readOutHalfSize;}

	 double crysHalfTrans() const       {return _crystalHalfTrans;}
	 double crystalHalfLength() const   {return _crystalDepth;}
         double wrapperThickness() const    {return _wrapperThickness; }
	 double shellThickness() const      {return _shellThickness;}

	 double getNonuniformity() const    {return _nonUniformity; }
	 double getTimeGap() const          {return _timeGap; }
	 double getElectronEdep() const     {return _electronEdep; }
	 double getElectronEmin() const     {return _electronEmin; }




         	 
 	 	 



       protected:

	  unsigned int _nDisk;
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
