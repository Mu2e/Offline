//
// Track direction convention 
//
// $Id: TrkFitDirection.cc,v 1.1 2012/07/25 20:56:57 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/25 20:56:57 $
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
