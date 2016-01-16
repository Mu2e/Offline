#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh
//
// $Id: HexMap.hh,v 1.6 2013/10/21 20:53:36 murat Exp $
// $Author: echenard $
// $Date: 2013/10/21 20:53:36 $
//

#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include <vector>
#include "CLHEP/Vector/TwoVector.h"



namespace mu2e {


    class VaneMapper : public CrystalMapper {


	public:

	    VaneMapper(int nCryX, int nCryY);

            int               nCrystalMax(int maxRing)                    {return _nCrystalX*_nCrystalY;}
            int               nApex()                               const {return _apexX.size();}
            double            apexX(int i)                          const {return _apexX.at(i);}
            double            apexY(int i)                          const {return _apexY.at(i);}

	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)            const;
            int               indexFromXY(double x, double y)       const;

	    std::vector<int>  neighbors(int thisIndex, unsigned int level=1) const;


	private:

	    std::vector<double>  _apexX;
	    std::vector<double>  _apexY;
	    int _nCrystalX;
	    int _nCrystalY;
	    
    };
}


#endif

