#ifndef CalorimeterGeom_SqaureShiftMapper_hh
#define CalorimeterGeom_SqaureShiftMapper_hh

#include "CalorimeterGeom/inc/CrystalMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include <vector>



namespace mu2e {

    class SquShiftLK {        

         public:

            SquShiftLK()             : l_(0),k_(0) {}
	    SquShiftLK(int l, int k) : l_(l),k_(k) {}

            void add(const SquShiftLK &x) {l_+=x.l_;k_+=x.k_;}

            int l_; 
	    int k_;                                
    };
    

    class SquareShiftMapper : public CrystalMapper {

	public:

	    SquareShiftMapper();
            virtual ~SquareShiftMapper() {};

            virtual int               nCrystalMax(int maxRing)        const {return 3*maxRing*(maxRing+1)+1;}
            virtual int               nApex()                         const {return apexX_.size();}
            virtual double            apexX(int i)                    const {return apexX_.at(i);}
            virtual double            apexY(int i)                    const {return apexY_.at(i);}

	    virtual CLHEP::Hep2Vector xyFromIndex(int thisIndex)          const;
            virtual int               indexFromXY(double x, double y)     const;
            virtual int               indexFromRowCol(int nRow, int nCol) const;

	    virtual std::vector<int>  neighbors(int thisIndex, unsigned int level=1) const;


	private:

	    SquShiftLK lk(int index)        const;
	    int index(const SquShiftLK& lk) const;
	    int ring(const SquShiftLK&lk)   const;

	    std::vector<SquShiftLK> step_;
	    std::vector<double>     apexX_;
	    std::vector<double>     apexY_;
    };
}

#endif

