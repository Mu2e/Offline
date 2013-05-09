#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// $Id: DiskCalorimeter.hh,v 1.5 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $
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
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {


    class DiskCalorimeter: public BaseCalorimeter{

       
       friend class DiskCalorimeterMaker;


       public:


          DiskCalorimeter()  {}
          ~DiskCalorimeter() {}

	  unsigned int nDisk()                const  {return _nSections;}
	  Disk const&  disk(int i)            const  {return static_cast<Disk const&>(section(i));}
	  double       diskSeparation(int i)  const  {return _diskSeparation.at(i);}

	  virtual double crystalHalfTrans()   const  {return _crystalHalfTrans;}
	  virtual double crystalHalfLength()  const  {return _crystalHalfLength;}
          virtual double crystalVolume()      const  {return 3.4641016*_crystalHalfTrans*_crystalHalfTrans*2.0*_crystalHalfLength;}


  	          bool             isInsideDisk(int idisk, CLHEP::Hep3Vector const& pos) const ;       	 
	  virtual bool             isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const ;       	 	 
	  virtual int              crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const ;
          virtual std::vector<int> neighbors(int crystalId, int level=1) const; 



        private:

	   std::vector<double> _diskInnerRadius;
	   std::vector<double> _diskOuterRadius;
	   std::vector<double> _diskSeparation;
	   std::vector<double> _diskRotAngle;  

	   double _crystalHalfTrans;
	   double _crystalHalfLength;

    };

}    

#endif /* CalorimeterGeom_DiskCalorimeter_hh */
