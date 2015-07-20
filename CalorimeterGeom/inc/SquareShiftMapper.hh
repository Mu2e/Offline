#ifndef CalorimeterGeom_SqaureShiftMapper_hh
#define CalorimeterGeom_SqaureShiftMapper_hh
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

    class SquShiftLK {        

         public:

            SquShiftLK(int l, int k) : _l(l),_k(k) {}

            int _l; 
	    int _k;                    

            void operator += (SquShiftLK x) {_l+=x._l;_k+=x._k;}
    };




    class SquareShiftMapper : public CrystalMapper {


	public:

	    SquareShiftMapper();

            int               nCrystalMax(int maxRing)        {return 3*maxRing*(maxRing+1)+1;}
            int               nApex()                   const {return _apexX.size();}
            double            apexX(int i)              const {return _apexX.at(i);}
            double            apexY(int i)              const {return _apexY.at(i);}

	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)      const;
            int               indexFromXY(double x, double y) const;

	    std::vector<int>  neighbors(int thisIndex, int level=1) const;



	private:

	    SquShiftLK   lk(int index)         const;
	    int          index(SquShiftLK& lk) const;
	    int          ring(SquShiftLK& lk)  const;

	    std::vector<SquShiftLK>   _step;
	    std::vector<double>       _apexX;
	    std::vector<double>       _apexY;
    };
}

#endif

