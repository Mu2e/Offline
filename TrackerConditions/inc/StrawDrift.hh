#ifndef TrackerConditions_StrawDrift_hh
#define TrackerConditions_StrawDrift_hh
#include <iostream>
#include <vector>
#include <string>
#include "TrackerConditions/inc/Types.hh"
//
// define the ends of a straw
namespace mu2e {
  struct point {
    float eField;
    float instVel;
  };

  class StrawDrift {
    public:
      StrawDrift (std::string filename);
      ~StrawDrift(){};

      std::vector<point> points;

      std::vector<float> distances;
      std::vector<float> effectiveVelocity;


  };
}
#endif

