#ifndef CalorimeterGeom_SquareMapper_hh
#define CalorimeterGeom_SquareMapper_hh
//
// Square position map generator:
//   tesselate a plane with squares starting from the center of the plane
//
//  original author : Bertrand Echenard (Caltech)
//
// Use basis vector, l and k, defined as
// l = right
// k = up
//
/*

 -----------------------
 |       |       |      |
 | -1 1  | 0 1   | 1 1  |
 |       |       |      |
 -----------------------
 |       |       |      |  l,k coordinates
 | -1 0  | 0 0   | 1 0  |
 |       |       |      |
 -----------------------
 |       |       |      |
 | -1 -1 | 0 -1  | 1 -1 |
 |       |       |      |
 -----------------------

  steps :  (1,0), (0,-1), (-1,0), (0,1) (clockwise from top left corner)
  segment: top=0, right=1, bottom=2,left=3

*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at -l,+l (top left corner),
//   then go n time each step to create the ring
//
// Neighbors add (+-1,0) or (0,+-1)
// next ring of neighbours, add +(-2,2) and go around the ring,...
//
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <array>


namespace mu2e {

  class SquLK {
    public:
      constexpr SquLK() = default;
      constexpr SquLK(int l, int k) : l_(l),k_(k) {}

      SquLK& operator+=(const SquLK& x) {l_ += x.l_; k_ += x.k_; return *this;}

      int l_;
      int k_;
  };


  class SquareMapper : public CrystalMapper {
    public:
      SquareMapper() = default;

      unsigned                  nCrystalMax(unsigned maxRing) const override;
      CLHEP::Hep2Vector         xyFromIndex(unsigned thisIndex) const override;
      unsigned                  indexFromXY(double x, double y) const override;
      unsigned                  indexFromRowCol(int row, int col) const override;
      int                       rowFromIndex(unsigned thisIndex) const override;
      int                       colFromIndex(unsigned thisIndex) const override;
      unsigned                  numNeighbors(unsigned level) const override;
      std::vector<unsigned>     neighbors(unsigned thisIndex,unsigned level) const override;


    private:
      SquLK    lk(unsigned index)     const;
      unsigned index(const SquLK& lk) const;
      unsigned ring(const  SquLK& lk) const;

      inline static constexpr std::array<SquLK,4>  step_{{{1,0},{0,-1},{-1,0},{0,1}}};
      inline static constexpr std::array<double,5> apexX_{{-0.5,0.5,0.5,-0.5,-0.5}};
      inline static constexpr std::array<double,5> apexY_{{-0.5,-0.5,0.5,0.5,-0.5}};
  };
}


#endif
