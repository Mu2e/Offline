#ifndef CalorimeterGeom_HexMapTest_hh
#define CalorimeterGeom_HexMapTest_hh


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

	  std::vector<int> getNeighbors(int index, int level=1) const;
	  CLHEP::Hep2Vector getXYPosition(int index) const;
          int getIndexFromXY(double x, double y) const;


	private:

	  std::vector<HexLK>  _step;

	  HexLK getLK(int index) const;
	  int getIndex(HexLK& lk) const;
	  int getRing(HexLK& lk) const;


    };

}
#endif

