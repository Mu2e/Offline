//
// reader to simulate a mock-up database
//
// Original author B Echenard 
//

#ifndef CalorimeterGeom_CrystalCondReader_hh
#define CalorimeterGeom_CrystalCondReader_hh

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>

namespace mu2e {

     class CrystalCondReader {

	  public:

             CrystalCondReader(const std::string& filename);
             ~CrystalCondReader() {};

             int                     nCrystal()               const;
             const CLHEP::Hep3Vector position(int index)      const; 
             const CLHEP::Hep3Vector size(int index)          const; 
             const void              print(std::ostream &os)  const;


	 private:
	 
	     void read(const std::string& filename);
             
             std::map<int, std::vector<double>> dataPosition_;
             std::map<int, std::vector<double>> dataSize_;

     };

}


#endif 
