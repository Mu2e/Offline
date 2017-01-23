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
	      origin_(),center_(),trackerCenter_(),diskInnerRadius_(),diskOuterRadius_(),
              diskSeparation_(),diskRotAngle_()   	      
	   {}
	     
           ~CaloGeomInfo() {}

           const CLHEP::Hep3Vector& origin()             const {return origin_;}
           const CLHEP::Hep3Vector& center()             const {return center_;}
           const CLHEP::Hep3Vector& trackerCenter()      const {return trackerCenter_;}
           
	   const std::vector<double>&  diskInnerRadius() const {return diskInnerRadius_;}
	   const std::vector<double>&  diskOuterRadius() const {return diskOuterRadius_;}
	   const std::vector<double>&  diskSeparation()  const {return diskSeparation_;}
	   const std::vector<double>&  diskRotAngle()    const {return diskRotAngle_;}



           void origin(const CLHEP::Hep3Vector vec)              {origin_ = vec;}
           void center(const CLHEP::Hep3Vector vec)              {center_ = vec;}
           void trackerCenter(const CLHEP::Hep3Vector vec)       {trackerCenter_ = vec;}

           void diskInnerRadius(const std::vector<double>&  vec) {diskInnerRadius_ = vec;}
           void diskOuterRadius(const std::vector<double>&  vec) {diskOuterRadius_ = vec;}
           void diskSeparation(const std::vector<double>&  vec)  {diskSeparation_ = vec;}
           void diskRotAngle(const std::vector<double>&  vec)    {diskRotAngle_ = vec;}




       private:

	  CLHEP::Hep3Vector origin_;
          CLHEP::Hep3Vector center_;
	  CLHEP::Hep3Vector trackerCenter_;

	  std::vector<double> diskInnerRadius_;
	  std::vector<double> diskOuterRadius_;
	  std::vector<double> diskSeparation_;
	  std::vector<double> diskRotAngle_;  


     };
}    

#endif 
