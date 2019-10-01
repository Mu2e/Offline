#ifndef ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
//
// $Id: EMFSimHitHistograms.hh,v 1.4 2013/07/30 18:45:00 wieschie Exp $
// $Author: wieschie $
// $Date: 2013/07/30 18:45:00 $
//
// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke
//

#include <string>
#include <map>

#include "boost/noncopyable.hpp"

#include "art_root_io/TFileDirectory.h"

#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

class TH1D;
class TH2D;

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
    void book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir);

    void fill(const ExtMonFNAL::ExtMon& extmon, const ExtMonFNALSimHitCollection& clusters);

  private:
    TH2D *hitTimes_;
    TH1D *energyDeposit_;
    TH1D *moduleHits_;
    std::map<ExtMonFNALModuleId, TH2D*> hitPosition_;
  };

} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh */
