//
// Track direction convention
//
//
// Original author David Brown, LBNL
//
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

namespace mu2e
{
  TrkFitDirection::TrkFitDirection(FitDirection fdir) : _fdir(fdir)
  {}

  std::string const&
  TrkFitDirection::name() const {
    switch (_fdir) {
      case upstream : {
        static std::string upname("Upstream");
        return upname;
      }
      case downstream: {
        static std::string downname("Downstream");
        return downname;
      }
    default : {
        static std::string unknown("Unknown");
        return unknown;
      }
    }
  }

  TrkFitDirection::FitDirection TrkFitDirection::fitDirectionFromName(std::string name) {
    TrkFitDirection::FitDirection fdir(TrkFitDirection::FitDirection::unknown);
    // convert to lowercase to protect against case-based issues
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c){ return std::tolower(c); });

    if     (name == "downstream") fdir = TrkFitDirection::FitDirection::downstream;
    else if(name == "upstream"  ) fdir = TrkFitDirection::FitDirection::upstream;
    return fdir;
  }

}
