#include "TrkChargeRecoOld/inc/PeakFitParams.hh"

namespace mu2e {
  namespace TrkChargeRecoOld {

    std::vector<std::string> PeakFitParams::_pnames { "earlyCharge",
						      "pedestal",
						      "time",
						      "charge",
						      "width",
						      "lateShift",
						      "lateCharge" };
  }
}
