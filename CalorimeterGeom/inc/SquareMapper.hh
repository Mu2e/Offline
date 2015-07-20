#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh
//
// $Id: HexMap.hh,v 1.6 2013/10/21 20:53:36 murat Exp $
// $Author: murat $
// $Date: 2013/10/21 20:53:36 $
//

//C++ // CLHEP includes
#include <vector>
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes
#include "CalorimeterGeom/inc/CrystalMapper.hh"


namespace mu2e {

    class SquLK {        

         public:

            SquLK(int l, int k) : _l(l),_k(k) {}

            int _l; 
	    int _k;                    

            void operator += (SquLK x) {_l+=x._l;_k+=x._k;}
    };




    class SquareMapper : public CrystalMapper {


	public:

	    SquareMapper();

            int               nCrystalMax(int maxRing)              {return (2*maxRing+1)*(2*maxRing+1);}
            int               nApex()                         const {return _apexX.size();}
            double            apexX(int i)                    const {return _apexX.at(i);}
            double            apexY(int i)                    const {return _apexY.at(i);}

	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)      const;
            int               indexFromXY(double x, double y) const;

	    std::vector<int>  neighbors(int thisIndex, int level=1) const;


	private:

	    SquLK    lk(int index)    const;
	    int      index(SquLK& lk) const;
	    int      ring(SquLK& lk)  const;

	    std::vector<SquLK>   _step;
	    std::vector<double>  _apexX;
	    std::vector<double>  _apexY;
    };
}


#endif

