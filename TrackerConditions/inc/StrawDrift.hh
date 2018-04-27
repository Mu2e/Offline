// Original author Jason Bono
//Feb 2018


#ifndef TrackerConditions_StrawDrift_hh
#define TrackerConditions_StrawDrift_hh
#include <iostream>
#include <vector>
#include <string>
#include "TrackerConditions/inc/Types.hh"


namespace mu2e {
  class StrawDrift {
  public:
    StrawDrift() : _initialized(false){};
    ~StrawDrift(){};
  
    void Initialize(std::string filename, float wirevoltage, int phiBins, int dIntegrationBins, float Bz);

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
    
    double GetAverageSpeed(double dist); //grabs the average nominal drift speed (phi = 0)
    double GetInstantSpeedFromT(double time); // (at phi = 0)
    double GetInstantSpeedFromD(double dist); // (at phi = 0)

    double GetGammaFromT(double time, double phi);
    double GetGammaFromD(double dist, double phi);
    double GetEffectiveSpeed(double dist, double phi); //grabs the lorentz corrected radial component of avg speed
    double D2T(double dist, double phi);
    double T2D(double time, double phi);

    bool isInitialized(){return _initialized;};
    
  private:
    double ConstrainAngle(double phi);

    bool _initialized;
    
    // 2-D array in distance and phi 
    std::vector<D2Tinfo> D2Tinfos;
    
    // 1-D array of speed with distance
    std::vector<float> distances;
    std::vector<float> instantSpeeds; // the instantaneous "nominal" speed
    std::vector<float> averageSpeeds; // the average "nominal" speed

    int _phiBins;
    
  };
}
#endif

