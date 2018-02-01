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
      //StrawDrift (std::string filename);
      StrawDrift (std::string filename, float wirevoltage); //JB
      ~StrawDrift(){};

      double integrateSpeed(float ddd); //this will build the table of average speed and distances
      double getEffectiveVelocity(double dist); //this will access the table of average speed and distances

      std::vector<point> points;

      
      std::vector<float> edistances; //JB
      std::vector<float> distances;
      std::vector<float> instantSpeeds; //JB
      std::vector<float> effectiveSpeeds; //JB
      //std::vector<float> effectiveVelocity;
      

  };
}
#endif

