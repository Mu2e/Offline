#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"

#include "cetlib_except/exception.h"

namespace mu2e {
  Mu2eG4TrajectoryControl::Mu2eG4TrajectoryControl(const Mu2eG4Config::TrajectoryControl_& tc)
    : produce_(tc.produce())
    , defaultMinPointDistance_{std::numeric_limits<double>::max() }
    , mcTrajectoryMinSteps_{std::numeric_limits<unsigned>::max() }
    , mcTrajectoryMomentumCut_{std::numeric_limits<double>::max() }
    , saveTrajectoryMomentumCut_{std::numeric_limits<double>::max() }
  {
    if(produce_) {

      if(!tc.defaultMinPointDistance(defaultMinPointDistance_)) {
        throw cet::exception("CONFIG")<< "Error: defaultMinPointDistance must be defined when MC Trajectory production is enabled\n";
      }

      if(!tc.mcTrajectoryMinSteps(mcTrajectoryMinSteps_)) {
        throw cet::exception("CONFIG")<< "Error: mcTrajectoryMinSteps must be defined when MC Trajectory production is enabled\n";
      }

      if(!tc.mcTrajectoryMomentumCut(mcTrajectoryMomentumCut_)) {
        throw cet::exception("CONFIG")<< "Error: mcTrajectoryMomentumCut must be defined when MC Trajectory production is enabled\n";
      }

      if(!tc.saveTrajectoryMomentumCut(saveTrajectoryMomentumCut_)) {
        throw cet::exception("CONFIG")<< "Error: saveTrajectoryMomentumCut must be defined when MC Trajectory production is enabled\n";
      }


      fhicl::ParameterSet volumeCutsPS;
      if(tc.perVolumeMinDistance.get_if_present(volumeCutsPS)) {
        const std::vector<std::string> volnames{volumeCutsPS.get_names()};
        for(const auto& k: volnames) {
          perVolumeMinDistance_[k] = volumeCutsPS.get<double>(k);
        }
      }
    }
  }
}
