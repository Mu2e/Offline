#ifndef CalorimeterGeom_SqaureShiftMapper_hh
#define CalorimeterGeom_SqaureShiftMapper_hh

#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include <vector>
#include "CLHEP/Vector/TwoVector.h"



namespace mu2e {

    class SquShiftLK {        

         public:

            SquShiftLK()             : _l(0),_k(0) {}
	    SquShiftLK(int l, int k) : _l(l),_k(k) {}

            void add(const SquShiftLK &x) {_l+=x._l;_k+=x._k;}

            int _l; 
	    int _k;                    
            
    };
    
    




    class SquareShiftMapper : public CrystalMapper {


	public:

	    SquareShiftMapper();

            int               nCrystalMax(int maxRing)        {return 3*maxRing*(maxRing+1)+1;}
            int               nApex()                   const {return _apexX.size();}
            double            apexX(int i)              const {return _apexX.at(i);}
            double            apexY(int i)              const {return _apexY.at(i);}

	    CLHEP::Hep2Vector xyFromIndex(int thisIndex)       const;
            int               indexFromXY(double x, double y)  const;

	    std::vector<int>  neighbors(int thisIndex, unsigned int level=1) const;



	private:

	    SquShiftLK   lk(int index)      const;
	    int index(SquShiftLK const &lk) const;
	    int ring(SquShiftLK const &lk)  const;

	    std::vector<SquShiftLK>   _step;
	    std::vector<double>       _apexX;
	    std::vector<double>       _apexY;
    };
}

#endif

