#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh


#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include <vector>
#include "CLHEP/Vector/TwoVector.h"



namespace mu2e {

    class SquLK {        

         public:

            SquLK()             : l_(0),k_(0) {}
            SquLK(int l, int k) : l_(l),k_(k) {}

            void add(const SquLK &x) {l_+=x.l_;k_+=x.k_;}

            int l_; 
	    int k_;                    
    };




    class SquareMapper : public CrystalMapper {

	public:

	    SquareMapper();
            virtual ~SquareMapper() {};

            virtual int               nCrystalMax(int maxRing)              const {return (2*maxRing+1)*(2*maxRing+1);}
            virtual int               nApex()                               const {return apexX_.size();}
            virtual double            apexX(int i)                          const {return apexX_.at(i);}
            virtual double            apexY(int i)                          const {return apexY_.at(i);}

	    virtual CLHEP::Hep2Vector xyFromIndex(int thisIndex)            const;
            virtual int               indexFromXY(double x, double y)       const;

	    virtual std::vector<int>  neighbors(int thisIndex, unsigned int level=1) const;


	private:

	    SquLK        lk(int index) const;
	    int index(const SquLK &lk) const;
	    int ring(const  SquLK &lk) const;

	    std::vector<SquLK>   step_;
	    std::vector<double>  apexX_;
	    std::vector<double>  apexY_;
    };
}


#endif

