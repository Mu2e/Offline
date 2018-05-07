//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////


#include "CalPatRec/inc/DeltaFinder2_types.hh"

namespace mu2e {
  namespace DeltaFinder2Types {
    
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result) {
      double x1, y1, x2, y2, nx1, ny1, nx2, ny2;
    
      const CLHEP::Hep3Vector& p1 = Hd1->fStraw->getMidPoint();

      x1 =  p1.x();
      y1 =  p1.y();

      const CLHEP::Hep3Vector& p2 = Hd2->fStraw->getMidPoint();
      x2 =  p2.x();
      y2 =  p2.y();

      const CLHEP::Hep3Vector& wdir1 = Hd1->fStraw->getDirection();
      nx1 = wdir1.x();
      ny1 = wdir1.y();

      const CLHEP::Hep3Vector& wdir2 = Hd2->fStraw->getDirection();
      nx2 = wdir2.x();
      ny2 = wdir2.y();

      double n1n2  = nx1*nx2+ny1*ny2;
      double r12n1 = (x1-x2)*nx1+(y1-y2)*ny1;
      double r12n2 = (x1-x2)*nx2+(y1-y2)*ny2;
      //-----------------------------------------------------------------------------
      // t1 and t2 are distances to the intersection point from the centers of the 
      // corresponding wires
      //-----------------------------------------------------------------------------
      Result->t1 = (n1n2*r12n2-r12n1)/(1-n1n2*n1n2);
      Result->t2 = (r12n2-n1n2*r12n1)/(1-n1n2*n1n2);

      // in 2D, the lines intersect, take one
      Result->x = x1+nx1*Result->t1;
      Result->y = y1+ny1*Result->t1;
      //-----------------------------------------------------------------------------
      // now define distances to the hits
      //-----------------------------------------------------------------------------
      const XYZVec* h1 = &Hd1->fPos->pos();
      Result->wd1 = (h1->x()-Result->x)*nx1+(h1->y()-Result->y)*ny1;
      const XYZVec* h2 = &Hd2->fPos->pos();
      Result->wd2 = (h2->x()-Result->x)*nx2+(h2->y()-Result->y)*ny2;

      return 0;
    }
  }
}
