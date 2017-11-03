#ifndef CalorimeterGeom_CaloGeomInfo_hh
#define CalorimeterGeom_CaloGeomInfo_hh
//
// Contains gometry info of disks
//
// Original author B. Echenard
//


#include "CLHEP/Vector/ThreeVector.h"
#include <vector>


namespace mu2e {

    class CaloGeomInfo {


       public:

           CaloGeomInfo(): 
	      origin_(),center_(),trackerCenter_(),diskCaseRadiusIn_(),diskCaseRadiusOut_(),
              diskZMotherShift_(),diskRotAngle_()   	      
	   {}
	     
           ~CaloGeomInfo() {}

           const CLHEP::Hep3Vector& origin()               const {return origin_;}
           const CLHEP::Hep3Vector& frontFace()            const {return frontFace_;}
           const CLHEP::Hep3Vector& trackerCenter()        const {return trackerCenter_;}
           
	   const std::vector<double>&  diskCaseRadiusIn()  const {return diskCaseRadiusIn_;}
	   const std::vector<double>&  diskCaseRadiusOut() const {return diskCaseRadiusOut_;}
	   const std::vector<double>&  diskZMotherShift()  const {return diskZMotherShift_;}
	   const std::vector<double>&  diskRotAngle()      const {return diskRotAngle_;}


           void origin(const CLHEP::Hep3Vector vec)                {origin_ = vec;}
           void frontFace(const CLHEP::Hep3Vector vec)             {frontFace_ = vec;}
           void trackerCenter(const CLHEP::Hep3Vector vec)         {trackerCenter_ = vec;}

           void diskCaseRadiusIn(const std::vector<double>&  vec)  {diskCaseRadiusIn_ = vec;}
           void diskCaseRadiusOut(const std::vector<double>&  vec) {diskCaseRadiusOut_ = vec;}
           void diskZMotherShift(const std::vector<double>&  vec)  {diskZMotherShift_ = vec;}
           void diskRotAngle(const std::vector<double>&  vec)      {diskRotAngle_ = vec;}



       private:

	  CLHEP::Hep3Vector origin_;
	  CLHEP::Hep3Vector frontFace_;
          CLHEP::Hep3Vector center_;
	  CLHEP::Hep3Vector trackerCenter_;

	  std::vector<double> diskCaseRadiusIn_;
	  std::vector<double> diskCaseRadiusOut_;
	  std::vector<double> diskZMotherShift_;
	  std::vector<double> diskRotAngle_;  
     };
}    

#endif 
