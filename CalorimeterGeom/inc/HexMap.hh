#ifndef CalorimeterGeom_HexMapTest_hh
#define CalorimeterGeom_HexMapTest_hh
//
// $Id: HexMap.hh,v 1.4 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $
//

//C++ includes
#include <vector>

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {

    class HexLK {        

         public:

            HexLK(int l, int k) : _l(l),_k(k) {}

            void operator += (HexLK x) {_l+=x._l;_k+=x._k;}
            int _l; 
	    int _k;                    

    };




    class HexMap {


	public:

	    HexMap();

	    CLHEP::Hep2Vector xyPosition(int thisIndex) const;
            int               indexFromXY(double x, double y) const;
	    std::vector<int>  neighbors(int thisIndex, int level=1) const;


	private:

	    HexLK lk(int index)    const;
	    int   index(HexLK& lk) const;
	    int   ring(HexLK& lk)  const;

	    std::vector<HexLK>  _step;
    };

}
#endif

