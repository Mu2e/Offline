#ifndef CalorimeterGeom_CrystalMapper_hh
#define CalorimeterGeom_CrystalMapper_hh
//
// Interface for classes describing the layout of the crystals in the disk
//
// Original author B. Echenard
//


#include "CLHEP/Vector/TwoVector.h"
#include <vector>


namespace mu2e {
    

    class CrystalMapper {


	public:

	   //no constructor for this interface
	   virtual ~CrystalMapper() {};

           virtual int                nCrystalMax(int maxRing) = 0;

           virtual int                nApex()                                      const = 0;
           virtual double             apexX(int i)                                 const = 0;
           virtual double             apexY(int i)                                 const = 0;

	   virtual CLHEP::Hep2Vector  xyFromIndex(int thisIndex)                   const = 0;
           virtual int                indexFromXY(double x, double y)              const = 0;

	   virtual std::vector<int>   neighbors(int thisIndex, unsigned int level) const = 0;

    };

}
#endif 
