#ifndef Mu2eUtilities_LineSegmentPCA_hh
#define Mu2eUtilities_LineSegmentPCA_hh
//========================================================================
// This class represents a LineSegment defined by two points (Hep2Vectors)
// representing the two ends of the line segment.
// The method IntersectResult determines if two LineSegments intersect.
//
// Original author Hans Wenzel
//
//========================================================================
// CLHEP includes
#include "CLHEP/Vector/TwoVector.h"

namespace mu2e {
  class LineSegmentPCA
  {
  private:
    CLHEP::Hep2Vector begin_;
    CLHEP::Hep2Vector end_;
    CLHEP::Hep2Vector intersection;
  public:
    LineSegmentPCA(const CLHEP::Hep2Vector& begin, const CLHEP::Hep2Vector& end)
      : begin_(begin), end_(end) {}
    enum IntersectResult{ PARALLEL, COINCIDENT, NOT_INTERSECTING, INTERSECTING };
    IntersectResult Intersect(const LineSegmentPCA& other_line, CLHEP::Hep2Vector& intersection);
  };
} // namespace mu2e

#endif /* Mu2eUtilities_LineSegmentPCA_hh */
