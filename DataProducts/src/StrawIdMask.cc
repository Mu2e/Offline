#include "DataProducts/inc/StrawIdMask.hh"

namespace mu2e {
  uint16_t StrawIdMask::fieldMask(StrawIdMask::field fval) {
    switch (fval) {
      case plane :
	return StrawId::_planemsk;
      case panel :
	return StrawId::_panelmsk;
      case straw :
	return StrawId::_strawmsk;
      case station :
	return StrawId::_stationmsk;
      case layer :
	return StrawId::_layermsk;
      default:
	return 0;
    }
    return 0;
  }
}
