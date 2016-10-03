// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Digitization_inc_PixelCharge_hh
#define ExtinctionMonitorFNAL_Digitization_inc_PixelCharge_hh

#include <ostream>
#include <map>
#include <queue>

#include "canvas/Persistency/Common/Ptr.h"

#include "DataProducts/inc/ExtMonFNALPixelId.hh"

namespace mu2e {

  class SimParticle;

  namespace ExtMonFNAL {

    struct PixelTimedChargeDeposit {
      double time;
      double charge;
      art::Ptr<SimParticle> particle;

      PixelTimedChargeDeposit(double t, double c, const art::Ptr<SimParticle>& p)
        : time(t), charge(c), particle(p)
      {}

      // We accumulated deposits in a priority_queue, and want earlier times
      // to come out first.  Thus the inverted less-than definition:
      bool operator<(const PixelTimedChargeDeposit& b) const {
        return b.time < this->time;
      }
    };

    std::ostream& operator<<(std::ostream& os, const PixelTimedChargeDeposit& dep) {
      return os<<"PTD("<<dep.time<<", "<<dep.charge<<", "<<dep.particle<<")";
    }

    // Queue ordered by time
    typedef std::priority_queue<PixelTimedChargeDeposit> PixelChargeHistory;

    typedef std::map<ExtMonFNALPixelId,PixelChargeHistory> PixelChargeCollection;

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Digitization_inc_PixelCharge_hh*/
