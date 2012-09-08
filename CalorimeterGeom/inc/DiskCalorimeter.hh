#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Disk Calorimeter. In order to insulate this class from
// knowledge of databases etc, this class can not know
// how to make itself.
//

// Look at Mu2eG4/inc/constructDiskCalorimeter.cc 
// for definition of geometry

// Original author B. Echenard
//

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"

//
// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"


namespace mu2e {

    class DiskCalorimeter: public Calorimeter{

       friend class DiskCalorimeterMaker;


       public:

	 DiskCalorimeter();
	 ~DiskCalorimeter() {}

	 unsigned int nDisks(void) const        {return _nDisks;}
	 Disk const& getDisk(int i) const       {return _disks.at(i);}
	 double getDiskSeparation(int i) const  {return _diskSeparation[i];}
	 double getDiskThickness(void) const    {return _diskThickness;}

	 unsigned int nRO(void) const;
         unsigned int nROPerCrystal(void) const {return _nROPerCrystal;}
	 double roHalfThickness(void) const     {return _readOutHalfThickness;}
	 double roHalfSize(void) const          {return _readOutHalfSize;}

	 double crysHexsize(void) const         {return _crystalHexsize;}
	 double crysDepth(void) const           {return _crystalDepth;}
         double wrapperThickness(void) const    {return _wrapperThickness; }
	 double shellThickness(void) const      {return _shellThickness;}

	 double getNonuniformity(void) const    {return _nonUniformity; }
	 double getTimeGap(void) const          {return _timeGap; }
	 double getElectronEdep(void) const     {return _electronEdep; }
	 double getElectronEmin(void) const     {return _electronEmin; }


	 CLHEP::Hep3Vector const& getOrigin(void) const  {return _origin;}


	 CLHEP::Hep3Vector toCrystalFrame(int roid, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector getCrystalOriginById(int id) const;
	 CLHEP::Hep3Vector toDiskFrame(int idisk, CLHEP::Hep3Vector const& pos) const ;
	 CLHEP::Hep3Vector fromDiskFrame(int idisk, CLHEP::Hep3Vector const& pos) const;   
         	 
         //conversion of crystal <-> readout id, 
	 //readout_id = crystal_id*nRoPerCrystal ... crystal_id*nRoPerCrystal + nRoPerCrystal-1		 
	 int getCrystalByRO(int roid) const   {return (roid/_nROPerCrystal);}
	 int getROBaseByCrystal(int id) const {return (id*_nROPerCrystal);}
	 int getDiskId(int id) const;
         int getCrystalIdInMap(int id) const;


	 
	 	 

       protected:

	  unsigned int _nDisks;
	  double _diskThickness;          
	  std::vector<double> _diskInnerRadius;
	  std::vector<double> _diskOuterRadius;
	  std::vector<double> _diskSeparation;
	  std::vector<double> _diskRotAngle;  //rotation w.r.t z-axi


	  unsigned int  _nROPerCrystal;
	  double _readOutHalfThickness;
	  double _readOutHalfSize;

	  double _crystalHexsize;
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

    

} //namespace mu2e

#endif /* CalorimeterGeom_DiskCalorimeter_hh */
