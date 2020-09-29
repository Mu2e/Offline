#include "DataProducts/inc/StrawIdMask.hh"
#include "cetlib_except/exception.h"

namespace mu2e {
  StrawIdMask::StrawIdMask(std::string const& tomask) {
    if(0 ==tomask.compare(0,5,"plane")){
      *this = StrawIdMask(plane);
    } else if(0 == tomask.compare(0,5,"panel")){
      *this = StrawIdMask(panel);
    } else if(0 == tomask.compare(0,11,"uniquepanel")){
      *this = StrawIdMask(uniquepanel);
    } else if(0 == tomask.compare(0,5,"straw")){
      *this = StrawIdMask(straw);
    } else if(0 == tomask.compare(0,11,"uniquestraw")){
      *this = StrawIdMask(uniquestraw);
    } else {
      throw cet::exception("CONFIG")
	<< "strawIdMask: supported values: 'plane', 'panel', 'uniquepanel', 'straw', 'uniquestraw'"
	<< "  Input was: " << tomask
	<< "\n";
    }
  }

  std::string StrawIdMask::levelName(Level level) {
    switch (level) {
      case plane:
	return std::string("plane");
      case panel:
	return std::string("panel");
      case uniquepanel:
	return std::string("uniquepanel");
      case straw:
	return std::string("straw");
      case uniquestraw:
	return std::string("uniquestraw");
      case unknown: default:
	return std::string("unknown");
    }
  }

  uint16_t StrawIdMask::levelMask(Level level) {
    switch (level) {
      case plane :
	return StrawId::_planemsk;
      case panel :
	return StrawId::_panelmsk;
      case straw :
	return StrawId::_strawmsk;
      case uniquepanel :
	return StrawId::_panelmsk | StrawId::_planemsk;
      case uniquestraw :
	return StrawId::_strawmsk | StrawId::_panelmsk | StrawId::_planemsk;
      default:
	return 0;
    }
    return 0;
  }
}
