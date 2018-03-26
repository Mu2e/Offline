// Original author Jason Bono
//Feb 2018


#ifndef TrackerConditions_StrawDrift_hh
#define TrackerConditions_StrawDrift_hh
#include <iostream>
#include <vector>
#include <string>
#include "TrackerConditions/inc/Types.hh"


namespace mu2e {
  struct point { //filled by the field/velocity file read in
    float eField;
    float instVel;
  };
  
  struct D2Tinfo{//the info here forms the basis of d2t and t2d
    float gamma;
    float effectiveSpeed;
    float time;
    float phi;
    float distance;
    float instantaneousSpeed; //needed later for error calculations
  };
  
  class StrawDrift {
  public:
    StrawDrift (std::string filename, float wirevoltage, float Bz);
    ~StrawDrift(){};
    
    double getAverageSpeed(double dist); //grabs the average nominal drift speed
    //double getGamma(double dist, double phi);
    double GetGammaFromT(double time, double phi);
    double GetGammaFromD(double dist, double phi);
    double GetInstantSpeedFromT(double time);
    double GetInstantSpeedFromD(double dist);

    
    
    double getEffectiveSpeed(double dist, double phi); //grabs the lorentz corrected radial component of avg speed
    double D2T(double dist, double phi);
    double T2D(double time, double phi);
    
    std::vector<point> points;
    std::vector<D2Tinfo> D2Tinfos;
    
    std::vector<float> edistances; //the set of distances that correspond to a given E-field
    std::vector<float> distances;
    std::vector<float> instantSpeeds; // the instantaneous "nominal" speed
    std::vector<float> averageSpeeds; // the average "nominal" speed
    double ConstrainAngle(double phi);
    
  };
}
#endif

