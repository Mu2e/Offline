//
// Sqaure position map generator:
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

#include "Offline/CalorimeterGeom/inc/SquareMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <map>
#include <cmath>


namespace mu2e {


      SquareMapper::SquareMapper() :
         step_(),
         apexX_({-0.5,0.5,0.5,-0.5,-0.5}),
         apexY_({-0.5,-0.5,0.5,0.5,-0.5})
      {
          step_.push_back( SquLK( 1, 0) );  //right
          step_.push_back( SquLK( 0,-1) );  //down
          step_.push_back( SquLK(-1, 0) );  //left
          step_.push_back( SquLK( 0, 1) );  //up
      }

      //--------------------------------------------------------------------------------
      int SquareMapper::nCrystalMax(int maxRing) const {return (2*maxRing+1)*(2*maxRing+1);}


      //--------------------------------------------------------------------------------
      CLHEP::Hep2Vector SquareMapper::xyFromIndex(int thisIndex) const
      {
          SquLK thisLK = lk(thisIndex);
          return CLHEP::Hep2Vector(thisLK.l_,thisLK.k_);
      }


      int SquareMapper::indexFromXY(double x0, double y0) const
      {
          int l = int(std::abs(x0)+0.5);
          int k = int(std::abs(y0)+0.5);
          if (x0<0) l *= -1;
          if (y0<0) k *= -1;

          SquLK lk(l,k);
          return index(lk);
      }


      //--------------------------------------------------------------------------------
      int SquareMapper::indexFromRowCol(int nRow, int nCol) const
      {
          SquLK lk(nCol,nRow);
          return index(lk);
      }

      int SquareMapper::rowFromIndex(int thisIndex) const
      {
          SquLK thisLK = lk(thisIndex);
          return thisLK.k_;
      }

      int SquareMapper::colFromIndex(int thisIndex) const
      {
          SquLK thisLK = lk(thisIndex);
          return thisLK.l_;
      }


      //--------------------------------------------------------------------------------
      bool SquareMapper::isInsideCrystal(double x, double y, const CLHEP::Hep3Vector& pos,
                                         const CLHEP::Hep3Vector& size) const
      {
          return (std::abs(x-pos.x()) < 0.5*size.x()) && (std::abs(y-pos.y()) < 0.5*size.y());
      }

      //--------------------------------------------------------------------------------
      std::vector<int> SquareMapper::neighbors(int thisIndex, int level)  const
      {
          if (level<1) return std::vector<int>{};

          std::vector<int> thisNeighbour;

          SquLK init = lk(thisIndex);
          SquLK lk(init.l_ - level, init.k_ + level);

          for (size_t i=0;i<step_.size();++i)
          {
              for (int iseg=0;iseg<2*level;++iseg)
              {
                 lk.add(step_[i]);
                 thisNeighbour.push_back( index(lk) );
              }
          }
          return thisNeighbour;
      }



      //--------------------------------------------------------------------------------
      SquLK SquareMapper::lk(int thisIndex) const
      {
          if (thisIndex==0) return SquLK(0,0);

          int nRing = int(0.5*sqrt(thisIndex) + 0.5);

          int nSeg  = (thisIndex - (2*nRing-1)*(2*nRing-1)) / (2*nRing);
          int nPos  = (thisIndex - (2*nRing-1)*(2*nRing-1)) % (2*nRing);

          if (nSeg==0) return SquLK( -nRing+nPos ,  nRing     );
          if (nSeg==1) return SquLK( nRing       ,  nRing-nPos);
          if (nSeg==2) return SquLK( nRing-nPos  , -nRing     );
          return SquLK(-nRing       , -nRing+nPos);
      }

      //--------------------------------------------------------------------------------
      int SquareMapper::index(const SquLK &thisLK) const
      {
          if (thisLK.l_==0 && thisLK.k_==0) return 0;

          int nRing = ring(thisLK);
          int pos   = (2*nRing-1)*(2*nRing-1);

          //add position along segment
          if ( thisLK.k_ ==  nRing && thisLK.l_ < nRing)   pos +=           nRing + thisLK.l_;
          if ( thisLK.l_ ==  nRing && thisLK.k_ > -nRing)  pos += 2*nRing + nRing - thisLK.k_;
          if ( thisLK.k_ == -nRing && thisLK.l_ > -nRing)  pos += 4*nRing + nRing - thisLK.l_;
          if ( thisLK.l_ == -nRing && thisLK.k_ < nRing)   pos += 6*nRing + nRing + thisLK.k_;
          return pos;
      }

      //--------------------------------------------------------------------------------
      int SquareMapper::ring(const SquLK &thisLK) const
      {
          return std::max(std::abs(thisLK.l_),std::abs(thisLK.k_));
      }

      //--------------------------------------------------------------------------------
      int SquareMapper::numNeighbors(int level) const
      {
          return 2*level*step_.size();
      }

}
