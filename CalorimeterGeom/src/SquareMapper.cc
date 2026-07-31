#include "Offline/CalorimeterGeom/inc/SquareMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <cmath>
#include <cstdlib>

namespace mu2e {

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::nCrystalMax(unsigned maxRing) const
  {
    return (2u*maxRing+1u)*(2u*maxRing+1u);
  }

  //--------------------------------------------------------------------------------
  CLHEP::Hep2Vector SquareMapper::xyFromIndex(unsigned thisIndex) const
  {
     SquLK thisLK = lk(thisIndex);
     return CLHEP::Hep2Vector(thisLK.l_,thisLK.k_);
  }

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::indexFromXY(double x0, double y0) const
  {
    int l = int(std::abs(x0)+0.5);
    int k = int(std::abs(y0)+0.5);

    if (x0 < 0) l = -l;
    if (y0 < 0) k = -k;

    return index({l,k});
  }

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::indexFromRowCol(int row, int col) const
  {
    return index({col,row});
  }

  //--------------------------------------------------------------------------------
  int SquareMapper::rowFromIndex(unsigned thisIndex) const
  {
    return lk(thisIndex).k_;
  }

  //--------------------------------------------------------------------------------
  int SquareMapper::colFromIndex(unsigned thisIndex) const
  {
    return lk(thisIndex).l_;
  }

  //--------------------------------------------------------------------------------
  std::vector<unsigned>
  SquareMapper::neighbors(unsigned thisIndex, unsigned level) const
  {
    if (!level)
      return {};

    std::vector<unsigned> result;
    result.reserve(numNeighbors(level));

    auto current = lk(thisIndex);
    current.l_ -= int(level);
    current.k_ += int(level);

    for (const auto& s : step_)
      for (unsigned j=0; j<2u*level; ++j) { current += s; result.push_back(index(current)); }
    return result;
  }

  //--------------------------------------------------------------------------------
  SquLK SquareMapper::lk(unsigned thisIndex) const
  {
    if (!thisIndex) return {0,0};

    const unsigned nRing = unsigned(0.5*sqrt(double(thisIndex))+0.5);
    const unsigned first = (2u*nRing-1u)*(2u*nRing-1u);
    const unsigned nSeg = (thisIndex-first)/(2u*nRing);
    const unsigned nPos = (thisIndex-first)%(2u*nRing);

    const int r = int(nRing);
    const int p = int(nPos);
    if (nSeg==0) return {-r+p,  r};
    if (nSeg==1) return { r,    r-p};
    if (nSeg==2) return { r-p, -r};
    return {-r,   -r+p};
  }

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::index(const SquLK& thisLK) const
  {
    if (thisLK.l_==0 && thisLK.k_==0) return 0u;

    unsigned nRing = ring(thisLK);
    unsigned pos = (2u*nRing-1u)*(2u*nRing-1u);
    int r = int(nRing);

    // add position along segment
    if (thisLK.k_ ==  r && thisLK.l_ <  r) pos += unsigned(r + thisLK.l_);
    if (thisLK.l_ ==  r && thisLK.k_ > -r) pos += unsigned(3*r - thisLK.k_);
    if (thisLK.k_ == -r && thisLK.l_ > -r) pos += unsigned(5*r - thisLK.l_);
    if (thisLK.l_ == -r && thisLK.k_ <  r) pos += unsigned(7*r + thisLK.k_);
    return pos;
  }

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::ring(const SquLK& thisLK) const
  {
    return unsigned(std::max(std::abs(thisLK.l_), std::abs(thisLK.k_)));
  }

  //--------------------------------------------------------------------------------
  unsigned SquareMapper::numNeighbors(unsigned level) const
  {
    return 2u*level*unsigned(step_.size());
  }

}
