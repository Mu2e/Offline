#ifndef CalorimeterGeom_HexMapTest_hh
#define CalorimeterGeom_HexMapTest_hh
//
// $Id: HexMap.hh,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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
           int _l; int _k;                    

    };




    class HexMap {


	public:

	  HexMap();

	  std::vector<int> neighbors(int index0, int level=1) const;
	  CLHEP::Hep2Vector xyPosition(int index0) const;
          int indexFromXY(double x, double y) const;


	private:

	  std::vector<HexLK>  _step;

	  HexLK lk(int index0) const;
	  int index(HexLK& lk) const;
	  int ring(HexLK& lk) const;


    };

}
#endif

