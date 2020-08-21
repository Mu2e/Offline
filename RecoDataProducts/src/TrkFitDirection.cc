//
// Track direction convention 
//
//
// Original author David Brown, LBNL
//
#include "RecoDataProducts/inc/TrkFitDirection.hh"

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
}
