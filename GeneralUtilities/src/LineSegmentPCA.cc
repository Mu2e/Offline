//========================================================================
// This class represents a LineSegment defined by two points (Hep2Vectors)
// representing the two ends of the line segment.
// The method IntersectResult determines if two LineSegments intersect.
//
// Original author Hans Wenzel
//
//=======================================================================
#include "GeneralUtilities/inc/LineSegmentPCA.hh"
using namespace std;

namespace mu2e {
    LineSegmentPCA::IntersectResult LineSegmentPCA::Intersect(const LineSegmentPCA& other_line, CLHEP::Hep2Vector& intersection)
    {
      float denom = ((other_line.end_.y() - other_line.begin_.y())*(end_.x() - begin_.x())) -
        ((other_line.end_.x() - other_line.begin_.x())*(end_.y() - begin_.y()));

      float nume_a = ((other_line.end_.x() - other_line.begin_.x())*(begin_.y() - other_line.begin_.y())) -
        ((other_line.end_.y() - other_line.begin_.y())*(begin_.x() - other_line.begin_.x()));

      float nume_b = ((end_.x() - begin_.x())*(begin_.y() - other_line.begin_.y())) -
        ((end_.y() - begin_.y())*(begin_.x() - other_line.begin_.x()));

      if(denom == 0.0f)
        {
          if(nume_a == 0.0f && nume_b == 0.0f)
            {
              return COINCIDENT;
            }
          return PARALLEL;
        }

      float ua = nume_a / denom;
      float ub = nume_b / denom;

      if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
        {
          // Get the intersection point.
          intersection =CLHEP::Hep2Vector(begin_.x() + ua*(end_.x() - begin_.x()),begin_.y() + ua*(end_.y() - begin_.y()));
          return INTERSECTING;
        }

      return NOT_INTERSECTING;
    }


} // end namespace mu2e
