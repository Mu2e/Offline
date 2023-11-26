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
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include <cmath>


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
      constexpr static const char* cxname = {"StrawDrift"};

      StrawDrift(double cc, size_t phiBins, double deltaD, std::vector<double> distances_dbins, std::vector<double> instantSpeed_dbins, std::vector<double> times_dbins,
          double deltaT, std::vector<double> distances_tbins, std::vector<double> times_tbins) :
        ProditionsEntity(cxname),
        _cc(cc),
        _phiBins(phiBins),   _deltaPhi(M_PI_2/static_cast<double>(_phiBins-1)),
        _deltaD(deltaD), _distances_dbins(distances_dbins),
        _instantSpeed_dbins(instantSpeed_dbins), _times_dbins(times_dbins),
        _deltaT(deltaT), _distances_tbins(distances_tbins), _times_tbins(times_tbins) {}

      virtual ~StrawDrift() = default;

      double GetCC() const { return _cc;};
      double GetAverageSpeed(double dist) const; // avg nom. drift speed (phi = 0)
      double GetInstantSpeedFromT(double time) const; // (at phi = 0)
      double GetInstantSpeedFromD(double dist) const; // (at phi = 0)
      double D2T(double dist, double phi) const;
      double T2D(double time, double phi, bool nonnegative=true) const;

      void print(std::ostream& os) const;

    private:

      // fold into first quadrant assuming the function
      // has x-z and y-z plane symmetry
      double foldPhi(double phi) const;
      void phiRange(double phi, size_t prange[2]) const;
      size_t timeIndex(double time) const;
      size_t distIndex(double dist) const;
      void indexRange(double time, double phi, size_t irange[2]) const;

      double _cc;
      size_t _phiBins;
      double _deltaPhi;

      double _deltaD;
      std::vector<double> _distances_dbins; // distances between points in the model
      std::vector<double> _instantSpeed_dbins; // the instantaneous "nominal" speed, 1D vs distance
      std::vector<double> _times_dbins; // 2d array vs distance and phi

      double _deltaT;
      std::vector<double> _distances_tbins; // 2d array vs time and phi
      std::vector<double> _times_tbins; // times between points for T2D


  };
}
#endif
