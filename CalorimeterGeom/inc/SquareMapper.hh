#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh


#include "CalorimeterGeom/inc/CrystalMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>



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

            virtual int               nCrystalMax(int maxRing)            const {return (2*maxRing+1)*(2*maxRing+1);}
	    virtual CLHEP::Hep2Vector xyFromIndex(int thisIndex)          const;
            virtual int               indexFromXY(double x, double y)     const;
            virtual int               indexFromRowCol(int nRow, int nCol) const;
            virtual bool              isInsideCrystal(double x, double y, 
                                                      const CLHEP::Hep3Vector& pos, 
                                                      const CLHEP::Hep3Vector& size) const; 

	    virtual std::vector<int>  neighbors(int thisIndex, int level=1) const;
            virtual const std::vector<double>& apexX() const {return apexY_;}
            virtual const std::vector<double>& apexY() const {return apexY_;}


	private:

	    SquLK lk(int index)        const;
	    int index(const SquLK &lk) const;
	    int ring(const  SquLK &lk) const;

	    std::vector<SquLK>   step_;
	    std::vector<double>  apexX_;
	    std::vector<double>  apexY_;
    };
}


#endif

