#ifndef ExtinctionMonitorFNAL_Analyses_EMFRecoClusterHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFRecoClusterHistograms_hh
//
//
// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke
//

#include <string>
#include <vector>

#include "boost/noncopyable.hpp"

#include "art_root_io/TFileDirectory.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

class TH1D;
class TH2D;

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }

  class EMFRecoClusterHistograms :  private boost::noncopyable {
  public:

    EMFRecoClusterHistograms();

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    // The default is to use the current module's TFileDirectory
    void book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath="");

    // Book histograms in the specified TFileDirectory.
    void book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir);

    void fill(const ExtMonFNALRecoClusterCollection& clusters);

  private:
    const ExtMonFNAL::ExtMon *extmon_; // non-owning

    TH2D* perPlaneClusterMultiplicity_;
    TH2D* perClusterHitMultiplicity_;
    TH2D* clusterClock_;

    TH2D* clusterWidth_;

    // per sensor
    std::vector<TH2D*> clusterPosition_;
  };

} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFRecoClusterHistograms_hh */
