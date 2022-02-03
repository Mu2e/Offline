#ifndef ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
//
//
// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke
//

#include <string>
#include <map>

#include "boost/noncopyable.hpp"

#include "Offline/DataProducts/inc/ExtMonFNALModuleId.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"

class TH1D;
class TH2D;

namespace art { class TFileDirectory; }

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }

  class EMFSimHitHistograms :  private boost::noncopyable {
  public:

    EMFSimHitHistograms() : hitTimes_(), energyDeposit_() {}

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    // The default is to use the current module's TFileDirectory
    void book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath="");

    // Book histograms in the specified TFileDirectory.
    void book(const ExtMonFNAL::ExtMon& extmon, const art::TFileDirectory& tfdir);

    void fill(const ExtMonFNAL::ExtMon& extmon, const ExtMonFNALSimHitCollection& clusters);

  private:
    TH2D *hitTimes_;
    TH1D *energyDeposit_;
    TH1D *moduleHits_;
    std::map<ExtMonFNALModuleId, TH2D*> hitPosition_;
  };

} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh */
