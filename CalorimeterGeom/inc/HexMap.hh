#ifndef CalorimeterGeom_HexMapTest_hh
#define CalorimeterGeom_HexMapTest_hh
//
// $Id: HexMap.hh,v 1.6 2013/10/21 20:53:36 murat Exp $
// $Author: murat $
// $Date: 2013/10/21 20:53:36 $
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

	    int               index(int l, int k) {HexLK lk(l,k); return index(lk);}
	    int               l(int index)        {return lk(index)._l;}
	    int               k(int index)        {return lk(index)._k;}

	    std::vector<int>  neighbors(int thisIndex, int level=1) const;

	    HexLK lk(int index)    const;
	    int   index(HexLK& lk) const;
	    int   ring(HexLK& lk)  const;

	private:

	    std::vector<HexLK>  _step;
    };

}
#endif

