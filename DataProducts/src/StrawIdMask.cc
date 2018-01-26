#include "DataProducts/inc/StrawIdMask.hh"

namespace mu2e {
  uint16_t StrawIdMask::fieldMask(StrawIdMask::field fval) {
    switch (fval) {
      case plane :
	return StrawId2::_planemsk;
      case panel :
	return StrawId2::_panelmsk;
      case straw :
	return StrawId2::_strawmsk;
      case station :
	return StrawId2::_stationmsk;
      case layer :
	return StrawId2::_layermsk;
      default:
	return 0;
    }
    return 0;
  }
}
