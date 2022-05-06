//
// Square position map generator:
//   tesselate a plane with squares, every row shifted horizontaly by 0.5 square size, starting from the center of the plane
//
//  original author : Bertrand Echenard (Caltech)
//
//
// Use basis vector, l and k, defined as
// l = up right
// k = down right
//
/*

       --------------------
       |         |        |
       |  0 -1   |   1 0  |
       |         |        |
       |         |        |
 ------------------------------
 |         |         |        |
 |  -1 -1  |   0 0   |  1 1   |   l,k coordinates
 |         |         |        |
 |         |         |        |
 ------------------------------
       |         |        |
       |  -1 0   |  0 1   |
       |         |        |
       |         |        |
       --------------------

  steps :  (1,1) (0,1) (-1,0) (-1,-1) (0,-1) (1,0) (clockwise from top left corner)

*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at 0,-l (top left corner),
//   then go n time each step to create the ring
//
// Neighbors add (0,-1) and go around the ring
// next ring of neighbours, add (0,-2) and go around the ring,...
//

#include "Offline/CalorimeterGeom/inc/SquareShiftMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <map>
#include <cmath>



namespace mu2e {


      SquareShiftMapper::SquareShiftMapper() :
        step_(),
        apexX_({-0.5,0.5,0.5,-0.5,-0.5}),
        apexY_({-0.5,-0.5,0.5,0.5,-0.5})
      {
          step_.push_back( SquShiftLK(  1,  1) ); //right
          step_.push_back( SquShiftLK(  0 , 1) ); //down right
          step_.push_back( SquShiftLK( -1,  0) ); //down left
          step_.push_back( SquShiftLK( -1, -1) ); //left
          step_.push_back( SquShiftLK(  0, -1) ); //up left
          step_.push_back( SquShiftLK(  1,  0) ); //up right
      }





      //--------------------------------------------------------------------------------
      CLHEP::Hep2Vector SquareShiftMapper::xyFromIndex(int thisIndex) const
      {
          SquShiftLK thisLK = lk(thisIndex);
          return CLHEP::Hep2Vector( (thisLK.l_+thisLK.k_)/2.0, (thisLK.l_-thisLK.k_) );
      }

      int SquareShiftMapper::indexFromXY(double x0, double y0) const
      {
          int l,k;
          int ny = (y0>0) ? int(std::abs(y0)+0.5) : -int(std::abs(y0)+0.5);

          if (ny%2==0)
          {
              int nx = (x0>0) ? int(std::abs(x0)+0.5) : -int(std::abs(x0)+0.5);
              l = nx+ny/2;
              k = nx-ny/2;
          } else {
              int nx = (x0>0) ? int(std::abs(x0)) : (-int(std::abs(x0))-1);
              l = nx + (ny+1)/2;
              k = nx - (ny-1)/2;
          }

          SquShiftLK lk(l,k);
          return index(lk);
      }



      //--------------------------------------------------------------------------------
      int SquareShiftMapper::indexFromRowCol(int nRow, int nCol) const
      {
          int k = nRow/2-nRow + nCol;
          int l = nRow/2 + nCol;

          SquShiftLK lk(l,k);
          return index(lk);
      }


      //--------------------------------------------------------------------------------
      bool SquareShiftMapper::isInsideCrystal(double x, double y, const CLHEP::Hep3Vector& pos,
                                              const CLHEP::Hep3Vector& size) const
      {
          return (std::abs(x-pos.x()) < 0.5*size.x()) && (std::abs(y-pos.y()) < 0.5*size.y());
      }



      //--------------------------------------------------------------------------------
      std::vector<int> SquareShiftMapper::neighbors(int thisIndex, int level)  const
      {
          std::vector<int> thisNeighbour;
          thisNeighbour.reserve(12);

          SquShiftLK init = lk(thisIndex);
          SquShiftLK lk(init.l_, init.k_ - level);

          for (size_t i=0;i<step_.size();++i)
          {
              for (int iseg=0;iseg<level;++iseg)
              {
                 lk.add(step_[i]);
                 thisNeighbour.push_back( index(lk) );
              }
          }
          return thisNeighbour;
      }


      //--------------------------------------------------------------------------------
      SquShiftLK SquareShiftMapper::lk(int thisIndex) const
      {
         if (thisIndex==0) return SquShiftLK(0,0);

         int nRing = int(0.5+sqrt(0.25+float(thisIndex-1)/3.0));
         int nSeg  = (thisIndex - 3*nRing*(nRing-1)-1) / nRing;
         int nPos  = (thisIndex - 3*nRing*(nRing-1)-1) % nRing;

         int l =          nPos*step_[nSeg].l_;
         int k = -nRing + nPos*step_[nSeg].k_;

         for (int i=0;i<nSeg;++i)
         {
            l += step_[i].l_*nRing;
            k += step_[i].k_*nRing;
         }

         return SquShiftLK(l,k);
      }


      //--------------------------------------------------------------------------------
      int SquareShiftMapper::index(const SquShiftLK& thisLK) const
      {
         if (thisLK.l_==0 && thisLK.k_==0) return 0;

         int nRing = ring(thisLK);
         int pos   = 3*nRing*(nRing-1)+1;

         //add position along segment -- a drawing helps...
         if (thisLK.k_ == -nRing && thisLK.l_==0) return pos;
         if (thisLK.l_ ==  nRing)                 return pos + nRing   + thisLK.k_;
         if (thisLK.l_ == -nRing)                 return pos + 4*nRing + std::abs(thisLK.k_);
         if (thisLK.k_ ==  nRing)                 return pos + 3*nRing - thisLK.l_;
         if (thisLK.k_ == -nRing)                 return pos + 6*nRing - std::abs(thisLK.l_);
         if (thisLK.l_ > thisLK.k_)               return pos + thisLK.l_;
         return pos + 3*nRing +std::abs(thisLK.l_);
      }


      //--------------------------------------------------------------------------------
      int SquareShiftMapper::ring(const SquShiftLK& thisLK) const
      {
          if (thisLK.l_*thisLK.k_>0) return std::max(std::abs(thisLK.l_),std::abs(thisLK.k_));
          return std::abs(thisLK.l_-thisLK.k_);
      }

}
