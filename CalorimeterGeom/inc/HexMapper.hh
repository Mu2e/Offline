#ifndef CalorimeterGeom_HexMapper_hh
#define CalorimeterGeom_HexMapper_hh
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

    class HexLK {        

         public:

            HexLK(int l, int k) : _l(l),_k(k) {}

            int _l; 
	    int _k;                    

            void operator += (HexLK x) {_l+=x._l;_k+=x._k;}
    };




    class HexMapper : public CrystalMapper {


	public:

	    HexMapper();

            int               nCrystalMax(int maxRing)        {return 1 + 3*maxRing*(maxRing-1);}
            int               nApex()                   const {return _apexX.size();}
            double            apexX(int i)              const {return _apexX.at(i);}
            double            apexY(int i)              const {return _apexY.at(i);}


	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)            const;
            int               indexFromXY(double x, double y)       const;

	    std::vector<int>  neighbors(int thisIndex, int level=1) const;


	private:

	    HexLK lk(int index)    const;
	    int   index(HexLK& lk) const;
	    int   ring(HexLK& lk)  const;

	    std::vector<HexLK>   _step;
	    std::vector<double>  _apexX;
	    std::vector<double>  _apexY;
    };

}
#endif

