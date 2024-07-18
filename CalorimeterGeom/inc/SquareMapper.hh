#ifndef CalorimeterGeom_SqaureMapper_hh
#define CalorimeterGeom_SqaureMapper_hh

#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
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

            int               nCrystalMax    (int maxRing)          const override;
            CLHEP::Hep2Vector xyFromIndex    (int thisIndex)        const override;
            int               indexFromXY    (double x, double y)   const override;
            int               indexFromRowCol(int nRow, int nCol)   const override;
            int               rowFromIndex   (int thisIndex)        const override;
            int               colFromIndex   (int thisIndex)        const override;
            bool              isInsideCrystal(double x, double y, const CLHEP::Hep3Vector& pos,
                                              const CLHEP::Hep3Vector& size) const override;

            std::vector<int>  neighbors(int thisIndex, int level=1) const override;
            const std::vector<double>& apexX() const override {return apexY_;}
            const std::vector<double>& apexY() const override {return apexY_;}


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

