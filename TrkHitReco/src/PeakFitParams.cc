#include "TrkHitReco/inc/PeakFitParams.hh"

namespace mu2e {
  namespace TrkHitReco {

    std::vector<std::string> PeakFitParams::_pnames { "earlyCharge",
						      "pedestal",
						      "time",
						      "charge",
						      "width",
						      "lateShift",
						      "lateCharge" };
  }
}
