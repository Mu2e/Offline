#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"



#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  Mu2eG4TrajectoryControl::Mu2eG4TrajectoryControl(const fhicl::ParameterSet& pset)
    : produce_(pset.get<bool>("produce"))
    , defaultMinPointDistance_{produce_ ? pset.get<double>("defaultMinPointDistance") : std::numeric_limits<double>::max() }
    , mcTrajectoryMinSteps_{produce_ ? pset.get<unsigned>("mcTrajectoryMinSteps") : std::numeric_limits<unsigned>::max() }
    , mcTrajectoryMomentumCut_{produce_ ? pset.get<double>("mcTrajectoryMomentumCut") : std::numeric_limits<double>::max() }
    , saveTrajectoryMomentumCut_{produce_ ? pset.get<double>("saveTrajectoryMomentumCut") : std::numeric_limits<double>::max() }
  {
    if(produce_) {
      const fhicl::ParameterSet& volumeCutsPS{pset.get<fhicl::ParameterSet>("perVolumeMinDistance")};
      const std::vector<std::string> volnames{volumeCutsPS.get_names()};
      for(const auto& k: volnames) {
        perVolumeMinDistance_[k] = volumeCutsPS.get<double>(k);
      }
    }
  }
}
