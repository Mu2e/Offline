// Andrei Gaponenko, 2015

#ifndef Mu2eG4_Mu2eG4TrajectoryControl_hh
#define Mu2eG4_Mu2eG4TrajectoryControl_hh

#include <string>
#include <map>

namespace fhicl { class ParameterSet; }

namespace mu2e {

  class Mu2eG4TrajectoryControl {
  public:
    typedef std::map<std::string, double> PerVolumeDistanceMap;

    explicit Mu2eG4TrajectoryControl(const fhicl::ParameterSet& pset);

    bool produce() const { return produce_; }
    double defaultMinPointDistance() const { return defaultMinPointDistance_; }
    unsigned mcTrajectoryMinSteps() const { return mcTrajectoryMinSteps_; }
    double mcTrajectoryMomentumCut() const { return mcTrajectoryMomentumCut_; }
    double saveTrajectoryMomentumCut() const { return saveTrajectoryMomentumCut_; }
    const PerVolumeDistanceMap& perVolumeMinDistance() const { return perVolumeMinDistance_; }

  private:
    bool produce_;
    double defaultMinPointDistance_;
    unsigned mcTrajectoryMinSteps_;
    double mcTrajectoryMomentumCut_;
    double saveTrajectoryMomentumCut_;
    PerVolumeDistanceMap perVolumeMinDistance_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4TrajectoryControl_hh */
