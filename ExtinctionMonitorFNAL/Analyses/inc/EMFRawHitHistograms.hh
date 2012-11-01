#ifndef ExtinctionMonitorFNAL_Analyses_EMFRawHitHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFRawHitHistograms_hh
//
// $Id: EMFRawHitHistograms.hh,v 1.2 2012/11/01 23:36:19 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:36:19 $
//
// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke
//

#include <string>
#include <map>

#include "boost/noncopyable.hpp"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "DataProducts/inc/ExtMonFNALChipId.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"

class TH1D;
class TH2D;

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }

  class EMFRawHitHistograms :  private boost::noncopyable {
  public:

    explicit EMFRawHitHistograms(const fhicl::ParameterSet& pset);

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    // The default is to use the current module's TFileDirectory
    void book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath="");

    // Book histograms in the specified TFileDirectory.
    void book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir);

    void fill(const ExtMonFNALRawHitCollection& clusters);

  private:
    bool foldTimeToMicrobunch_;
    int numMicrobunchTicks_;

    TH1D* hitClock_;
    TH1D* hitToT_;

    // per-chip histograms
    std::map<ExtMonFNALChipId, TH2D*> chipOccupancy_;
  };

} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFRawHitHistograms_hh */
