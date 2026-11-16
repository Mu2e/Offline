#include "Offline/CalorimeterGeom/inc/SquareShiftMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <cmath>
#include <cstdlib>

namespace mu2e {

  //--------------------------------------------------------------------------------
  unsigned SquareShiftMapper::nCrystalMax(unsigned maxRing) const
  {
    return 3u*maxRing*(maxRing+1u)+1u;
  }

  //--------------------------------------------------------------------------------
  CLHEP::Hep2Vector SquareShiftMapper::xyFromIndex(unsigned thisIndex) const
  {
     SquShiftLK thisLK = lk(thisIndex);
     return CLHEP::Hep2Vector( (thisLK.l_+thisLK.k_)/2.0, (thisLK.l_-thisLK.k_) );
  }

  //--------------------------------------------------------------------------------
  unsigned SquareShiftMapper::indexFromXY(double x0, double y0) const
  {
    int l=0, k=0;
    int ny = (y0>0) ? int(std::abs(y0)+0.5) : -int(std::abs(y0)+0.5);

    if (ny%2==0) {
      int nx = (x0>0) ? int(std::abs(x0)+0.5) : -int(std::abs(x0)+0.5);
      l = nx + ny/2;
      k = nx - ny/2;
    } else {
      int nx = (x0>0) ? int(std::abs(x0)) : -int(std::abs(x0))-1;
      l = nx + (ny+1)/2;
      k = nx - (ny-1)/2;
    }

    return index({l,k});
  }

  unsigned SquareShiftMapper::indexFromRowCol(int row, int col) const
  {
    const int p = ((row % 2) + 2) % 2; // 0 if row even, 1 if odd (negatives safe)
    const int s = 2*col + p;           // s = l + k
    const int l = (s + row)/2;
    const int k = (s - row)/2;
    return index({l,k});
  }

  //--------------------------------------------------------------------------------
  int SquareShiftMapper::rowFromIndex(unsigned thisIndex) const
  {
    auto lk0 = lk(thisIndex);
    return lk0.l_-lk0.k_;
  }

  //--------------------------------------------------------------------------------
  int SquareShiftMapper::colFromIndex(unsigned thisIndex) const
  {
    auto lk0 = lk(thisIndex);
    int s = lk0.l_ + lk0.k_;
    return (s >= 0) ? s/2 : -(((-s) + 1)/2); // floor((l+k)/2)
  }

  //--------------------------------------------------------------------------------
  std::vector<unsigned> SquareShiftMapper::neighbors(unsigned thisIndex, unsigned level) const
  {
    if (!level) return {};

    std::vector<unsigned> result;
    result.reserve(numNeighbors(level));

    auto current = lk(thisIndex);
    current.k_ -= int(level);
    for (const auto& s : step_)
      for (unsigned j=0; j<level; ++j) { current += s; result.push_back(index(current)); }
    return result;
  }

  //--------------------------------------------------------------------------------
  SquShiftLK SquareShiftMapper::lk(unsigned thisIndex) const
  {
    if (!thisIndex) return {0,0};

    unsigned nRing = unsigned(0.5+std::sqrt(0.25+(thisIndex-1u)/3.0));
    unsigned nSeg = (thisIndex-3u*nRing*(nRing-1u)-1u)/nRing;
    unsigned nPos = (thisIndex-3u*nRing*(nRing-1u)-1u)%nRing;

    if (nSeg>=step_.size()) return {0,0};

    int l = int(nPos)*step_[nSeg].l_;
    int k = -int(nRing)+int(nPos)*step_[nSeg].k_;

    for (unsigned i=0; i<nSeg; ++i) {
      l += int(nRing)*step_[i].l_;
      k += int(nRing)*step_[i].k_;
    }

    return {l,k};
  }

  //--------------------------------------------------------------------------------
  unsigned SquareShiftMapper::index(const SquShiftLK& lk0) const
  {
    if (!lk0.l_ && !lk0.k_) return 0u;

    unsigned nRing = ring(lk0);
    unsigned pos = 3u*nRing*(nRing-1u)+1u;
    int r = int(nRing);

    if (lk0.k_==-r && lk0.l_==0) return pos;
    if (lk0.l_== r)    return pos + unsigned(r+lk0.k_);
    if (lk0.l_==-r)    return pos + unsigned(4*r+std::abs(lk0.k_));
    if (lk0.k_== r)    return pos + unsigned(3*r-lk0.l_);
    if (lk0.k_==-r)    return pos + unsigned(6*r-std::abs(lk0.l_));
    if (lk0.l_>lk0.k_) return pos + unsigned(lk0.l_);
    return pos + 3u*nRing + unsigned(std::abs(lk0.l_));
  }

  //--------------------------------------------------------------------------------
  unsigned SquareShiftMapper::ring(const SquShiftLK& lk0) const
  {
    if (lk0.l_*lk0.k_>0)
      return unsigned(std::max(std::abs(lk0.l_),std::abs(lk0.k_)));
    return unsigned(std::abs(lk0.l_-lk0.k_));
  }

  //--------------------------------------------------------------------------------
  unsigned SquareShiftMapper::numNeighbors(unsigned level) const
  {
    return level*unsigned(step_.size());
  }

}
