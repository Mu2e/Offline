#ifndef ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh
//
// $Id: EMFSimHitHistograms.hh,v 1.3 2012/11/01 23:43:17 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:43:17 $
//
// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke
//

#include <string>
#include <map>

#include "boost/noncopyable.hpp"

#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "DataProducts/inc/ExtMonFNALSensorId.hh"

class TH1D;
class TH2D;

namespace mu2e {

  class ExtMonFNALSimHitCollection;
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

    void fill(const ExtMonFNALSimHitCollection& clusters);

  private:
    TH2D *hitTimes_;
    TH1D *energyDeposit_;
    std::map<ExtMonFNALSensorId, TH2D*> hitPosition_;
  };

} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFSimHitHistograms_hh */
