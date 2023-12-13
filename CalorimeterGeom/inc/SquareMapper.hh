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

            virtual int               nCrystalMax    (int maxRing)          const override;
            virtual CLHEP::Hep2Vector xyFromIndex    (int thisIndex)        const override;
            virtual int               indexFromXY    (double x, double y)   const override;
            virtual int               indexFromRowCol(int nRow, int nCol)   const override;
            virtual int               rowFromIndex   (int thisIndex)        const override;
            virtual int               colFromIndex   (int thisIndex)        const override;
            virtual bool              isInsideCrystal(double x, double y,
                                                      const CLHEP::Hep3Vector& pos,
                                                      const CLHEP::Hep3Vector& size) const override;

            virtual std::vector<int>  neighbors(int thisIndex, int level=1) const override;
            virtual const std::vector<double>& apexX() const override {return apexY_;}
            virtual const std::vector<double>& apexY() const override {return apexY_;}


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

