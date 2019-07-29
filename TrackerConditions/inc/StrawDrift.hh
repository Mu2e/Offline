#ifndef TrackerConditions_StrawDrift_hh
#define TrackerConditions_StrawDrift_hh

//
// A conditions entity to hold a model of straw drift 
// (velocity as a function of field or distance)
// and provide it though various accessors
//


#include <iostream>
#include <vector>
#include <string>
#include "DataProducts/inc/TrkTypes.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"


namespace mu2e {
  class StrawDrift : public ProditionsEntity {
  public:

    struct D2Tinfo{//the info here forms the basis of d2t and t2d
      float gamma;
      float effectiveSpeed;
      float time;
      float phi;
      float distance;
      float instantaneousSpeed; //needed later for error calculations
    };

    typedef std::shared_ptr<StrawDrift> ptr_t;
    typedef std::shared_ptr<const StrawDrift> cptr_t;

    StrawDrift():_name("StrawDrift") {}
    StrawDrift( std::vector<D2Tinfo> D2Tinfos, std::vector<float> distances,
                std::vector<float> instantSpeeds, std::vector<float> averageSpeeds,
		int phiBins ) : _name("StrawDrift"),
      _D2Tinfos(D2Tinfos),  _distances(distances), 
      _instantSpeeds(instantSpeeds), _averageSpeeds(averageSpeeds),
      _phiBins(phiBins) {}

    virtual ~StrawDrift() {}

    double GetAverageSpeed(double dist) const; // avg nom. drift speed (phi = 0)
    double GetInstantSpeedFromT(double time) const; // (at phi = 0)
    double GetInstantSpeedFromD(double dist) const; // (at phi = 0)
    double GetGammaFromT(double time, double phi) const;
    double GetGammaFromD(double dist, double phi) const;
    // the lorentz corrected radial component of avg speed
    double GetEffectiveSpeed(double dist, double phi) const; 
    double D2T(double dist, double phi) const;
    double T2D(double time, double phi) const;

    void print(std::ostream& os) const;
    std::string const& name() const { return _name; }

  private:
    std::string _name;

    // fold into first quadrant assuming the function
    // has x-z and y-z plane symmetry
    double ConstrainAngle(double phi) const;
    // find distance bin
    size_t lowerDistanceBin(double dist) const;

    // 2-D array in distance and phi 
    std::vector<D2Tinfo> _D2Tinfos;
    
    // 1-D array of speed with distance
    std::vector<float> _distances; // distances between points in the model
    std::vector<float> _instantSpeeds; // the instantaneous "nominal" speed
    std::vector<float> _averageSpeeds; // the average "nominal" speed
    
    size_t _phiBins;
    
  };
}
#endif

