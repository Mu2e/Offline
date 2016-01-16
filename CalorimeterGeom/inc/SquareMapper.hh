#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh


#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include <vector>
#include "CLHEP/Vector/TwoVector.h"



namespace mu2e {

    class SquLK {        

         public:

            SquLK()             : _l(0),_k(0) {}
            SquLK(int l, int k) : _l(l),_k(k) {}

            void add(const SquLK &x) {_l+=x._l;_k+=x._k;}

            int _l; 
	    int _k;                    
    };




    class SquareMapper : public CrystalMapper {


	public:

	    SquareMapper();

            int               nCrystalMax(int maxRing)                    {return (2*maxRing+1)*(2*maxRing+1);}
            int               nApex()                               const {return _apexX.size();}
            double            apexX(int i)                          const {return _apexX.at(i);}
            double            apexY(int i)                          const {return _apexY.at(i);}

	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)            const;
            int               indexFromXY(double x, double y)       const;

	    std::vector<int>  neighbors(int thisIndex, unsigned int level=1) const;


	private:

	    SquLK        lk(int index) const;
	    int index(const SquLK &lk) const;
	    int ring(const  SquLK &lk) const;

	    std::vector<SquLK>   _step;
	    std::vector<double>  _apexX;
	    std::vector<double>  _apexY;
    };
}


#endif

