#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// $Id: DiskCalorimeter.hh,v 1.7 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
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

	  
	  
	  //disk components
	  unsigned int nDisk()                  const  {return _nSections;}
	  Disk const&  disk(int i)              const  {return static_cast<Disk const&>(section(i));}
	  double       diskSeparation(int i)    const  {return _diskSeparation.at(i);}

	  
	  //geometry components
	  virtual bool              isInsideCalorimeter(CLHEP::Hep3Vector const& pos)               const ;       	 	 
          virtual bool              isInsideSection(int iSection, CLHEP::Hep3Vector const& pos)     const ;

	  
	  //crystal id and neighbors component
	  virtual int               crystalIdxFromPosition(CLHEP::Hep3Vector const& pos)            const ;
          virtual std::vector<int>  neighborsByLevel(int crystalId, int level)                      const; 
	  virtual void              print()                                                         const;



        private:

	   std::vector<double> _diskInnerRadius;
	   std::vector<double> _diskOuterRadius;
	   std::vector<double> _diskSeparation;
	   std::vector<double> _diskRotAngle;  

    };

}    

#endif 
