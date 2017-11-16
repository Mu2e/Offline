#include "TrkChargeReco/inc/PeakFitParams.hh"

namespace mu2e {
  namespace TrkChargeReco {

    std::vector<std::string> PeakFitParams::_pnames { "earlyCharge",
						      "pedestal",
						      "time",
						      "charge",
						      "width",
						      "lateShift",
						      "lateCharge" };
  }
}
