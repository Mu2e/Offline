#ifndef CalorimeterGeom_SquareShiftMapper_hh
#define CalorimeterGeom_SquareShiftMapper_hh
//
// Square position map generator:
//   tesselate a plane with squares, every row shifted horizontaly by 0.5 square size,
//   starting from the center of the plane
//
//  original author : Bertrand Echenard (Caltech)
//
// Use basis vector, l and k, defined as
// l = up right
// k = down right
//
//       --------------------
//       |         |        |
//       |  0 -1   |   1 0  |
//       |         |        |
//       |         |        |
// ------------------------------
// |         |         |        |
// |  -1 -1  |   0 0   |  1 1   |   l,k coordinates
// |         |         |        |
// |         |         |        |
// ------------------------------
//       |         |        |
//       |  -1 0   |  0 1   |
//       |         |        |
//       |         |        |
//       --------------------
//
//  steps :  (1,1) (0,1) (-1,0) (-1,-1) (0,-1) (1,0) (clockwise from top left corner)
//
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at 0,-l (top left corner),
//   then go n time each step to create the ring
//
// Neighbors add (0,-1) and go around the ring
// next ring of neighbours, add (0,-2) and go around the ring,...
//
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <array>


namespace mu2e {

  class SquShiftLK {
    public:
      constexpr SquShiftLK() = default;
      constexpr SquShiftLK(int l, int k) : l_(l),k_(k) {}

      SquShiftLK& operator+=(const SquShiftLK& x) {l_ += x.l_; k_ += x.k_; return *this;}

      int l_;
      int k_;
  };


  class SquareShiftMapper : public CrystalMapper {
    public:
      SquareShiftMapper() = default;

      unsigned                  nCrystalMax(unsigned maxRing) const override;
      CLHEP::Hep2Vector         xyFromIndex(unsigned thisIndex) const override;
      unsigned                  indexFromXY(double x, double y) const override;
      unsigned                  indexFromRowCol(int row, int col) const override;
      int                       rowFromIndex(unsigned thisIndex) const override;
      int                       colFromIndex(unsigned thisIndex) const override;
      unsigned                  numNeighbors(unsigned level) const override;
      std::vector<unsigned>     neighbors(unsigned thisIndex,unsigned level) const override;


    private:
      SquShiftLK lk(unsigned index) const;
      unsigned index(const SquShiftLK& lk) const;
      unsigned ring(const SquShiftLK& lk) const;

      inline static constexpr std::array<SquShiftLK,6> step_{{{1,1},{0,1},{-1,0},{-1,-1},{0,-1},{1,0}}};
      inline static constexpr std::array<double,5> apexX_{{-0.5,0.5,0.5,-0.5,-0.5}};
      inline static constexpr std::array<double,5> apexY_{{-0.5,-0.5,0.5,0.5,-0.5}};
  };
}

#endif
