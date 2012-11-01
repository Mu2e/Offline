// Andrei Gaponenko, following GeneratorSummaryHistograms by Rob Kutschke

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFPatRecEffHistograms.hh"

#include <sstream>

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    EMFPatRecEffHistograms::EMFPatRecEffHistograms(unsigned cutMinCommonClusters)
      : cutMinCommonClusters_(cutMinCommonClusters)
      , extmon_()
      , hcommon_()
      , effMultiplicity_()
    {}

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void EMFPatRecEffHistograms::book(const ExtMon& extmon, const std::string& relativePath)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
      book (extmon, tfdir);
    }

    // Book the histograms.
    void EMFPatRecEffHistograms::book(const ExtMon& extmon, art::TFileDirectory& tfdir) {

      hcommon_ = tfdir.make<TH1D>("ncommon", "Num common clusters for best match", 7, -0.5, 6.5);

      effMultiplicity_ = tfdir.make<TEfficiency>("effMultiplicity", "Eff vs multiplicity;signal multiplicity;#epsilon",
                                                 40, 0.5, 200.5);

      effMultiplicity_->SetStatisticOption(TEfficiency::kFNormal);

    } // end EMFPatRecEffHistograms::book()

    //================================================================
    void EMFPatRecEffHistograms::Fillable::fill(unsigned denominatorParticleIndex) {

      // Find all matches
      const std::vector<const ExtMonFNALTrkMatchInfo*>& matchInfo = trackFinder_.data(denominatorParticleIndex);

      if(!matchInfo.empty()) {
        // select best match
        unsigned ibest = 0;
        for(unsigned i=1; i<matchInfo.size(); ++i) {
          if(matchInfo[ibest]->nCommonClusters() < matchInfo[i]->nCommonClusters()) {
            ibest = i;
          }
        }
        const ExtMonFNALTrkMatchInfo& bestMatch = *matchInfo[ibest];
        parent_->hcommon_->Fill(bestMatch.nCommonClusters());

        bool recoOK = (parent_->cutMinCommonClusters_ <= bestMatch.nCommonClusters() );
        parent_->effMultiplicity_->Fill(recoOK, multiplicity_);

      }
      else {
        parent_->hcommon_->Fill(0);
        parent_->effMultiplicity_->Fill(false, multiplicity_);
      }

    } // fill()

  } // end namespace ExtMonFNAL
} // end namespace mu2e
