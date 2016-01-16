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
//
// Look at Mu2eG4/inc/constructDiskCalorimeter.cc for definition of geometry
//
// Original author B. Echenard


#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>



namespace mu2e {


    class DiskCalorimeter: public BaseCalorimeter {

       
       friend class DiskCalorimeterMaker;

       public:


          DiskCalorimeter()  {}
          ~DiskCalorimeter() {}

	  
	  
	  //disk components
	  unsigned int nDisk()                  const  {return _nSections;}
	  Disk const&  disk(int i)              const  {return static_cast<Disk const&>(section(i));}
	  double       diskSeparation(int i)    const  {return _diskSeparation.at(i);}

	  
	  //geometry components
	  virtual bool              isInsideCalorimeter(const CLHEP::Hep3Vector &pos)                                 const;       	 	 
          virtual bool              isInsideSection(int iSection, const CLHEP::Hep3Vector &pos)                       const;
	  virtual bool              isContainedSection(const CLHEP::Hep3Vector &front, const CLHEP::Hep3Vector &back) const;

	  
	  //crystal id and neighbors component
	  virtual int               crystalIdxFromPosition(const CLHEP::Hep3Vector &pos)                    const;
          virtual std::vector<int>  neighborsByLevel(int crystalId, int level)                              const; 
	  virtual void              print(std::ostream &os = std::cout)                                     const;



        private:

	   std::vector<double> _diskInnerRadius;
	   std::vector<double> _diskOuterRadius;
	   std::vector<double> _diskSeparation;
	   std::vector<double> _diskRotAngle;  

    };

}    

#endif 
