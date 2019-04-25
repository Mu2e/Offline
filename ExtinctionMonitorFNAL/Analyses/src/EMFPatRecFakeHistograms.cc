// ExtMonFNAL fake rate analysis.
//
// Andrei Gaponenko

#include "ExtinctionMonitorFNAL/Analyses/inc/EMFPatRecFakeHistograms.hh"

#include <sstream>

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    EMFPatRecFakeHistograms::EMFPatRecFakeHistograms(unsigned cutMinCommonClusters)
      : cutMinCommonClusters_(cutMinCommonClusters)
      , extmon_()

      , fracDups_()
      , fracInteractions_()
      , fracFakes_()
    {}

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void EMFPatRecFakeHistograms::book(const ExtMon& extmon, const std::string& relativePath)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
      book (extmon, tfdir);
    }

    // Book the histograms.
    void EMFPatRecFakeHistograms::book(const ExtMon& extmon, art::TFileDirectory& tfdir) {

      fracDups_ = tfdir.make<TEfficiency>("fracDups", "Fraction of duplicates vs multiplicity;signal multiplicity",
                                          500, 0.5, 500.5);

      fracDups_->SetStatisticOption(TEfficiency::kFNormal);

//FIXME:      fracInteractions_ = tfdir.make<TEfficiency>("fracInteractions", "Fraction of interacted tracks vs multiplicity;signal multiplicity",
//FIXME:                                                 40, 0.5, 200.5);
//FIXME:
//FIXME:      fracInteractions_->SetStatisticOption(TEfficiency::kFNormal);

      fracFakes_ = tfdir.make<TEfficiency>("fracFakes", "Fraction of fake tracks vs multiplicity;signal multiplicity",
                                           500, 0.5, 500.5);

      fracFakes_->SetStatisticOption(TEfficiency::kFNormal);

    } // end EMFPatRecFakeHistograms::book()

    //================================================================
    int EMFPatRecFakeHistograms::Fillable::fill(unsigned denominatorTrackIndex) {

      // Find all matches
      const std::vector<const SimParticle*>& particles = particleFinder_.at(denominatorTrackIndex);
      const std::vector<const ExtMonFNALTrkMatchInfo*>& matchInfo = particleFinder_.data(denominatorTrackIndex);

      bool goodMatchFound = false;
      bool duplicateTrack = false;

      for(unsigned i=0; i<matchInfo.size(); ++i) {
        if(parent_->cutMinCommonClusters_ <= matchInfo[i]->nCommonClusters()) {
          goodMatchFound = true;
          if(!usedParticles_.insert(particles[i]).second) {
            duplicateTrack = true;
            goodMatchFound = false;
          }
        }
      }

      parent_->fracFakes_->Fill(!goodMatchFound, multiplicity_);
      parent_->fracDups_->Fill(duplicateTrack, multiplicity_);

      return !goodMatchFound;

    } // fill()

  } // end namespace ExtMonFNAL
} // end namespace mu2e
