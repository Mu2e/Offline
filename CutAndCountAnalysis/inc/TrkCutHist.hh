// Histograms filled during track selection: neither "all" nor
// "accepted" tracks
//
// Andrei Gaponenko, 2016

#ifndef CutAndCountAnalysis_inc_TrkCutHist_hh
#define CutAndCountAnalysis_inc_TrkCutHist_hh

#include <string>

#include "boost/noncopyable.hpp"

#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {
  namespace CutAndCount {

    struct TrkCutHist: private boost::noncopyable {
      explicit TrkCutHist(art::TFileDirectory tfdir, const std::string& relpath);
      TH1 *trkqual;
      TH1 *td;
      TH1 *d0;
      TH1 *rmax;
      TH1 *t0;
      TH1 *caloMatchChi2;
      TH1 *caloClusterEnergy;
      TH1 *momentum;
    private:
      art::TFileDirectory getdir(art::TFileDirectory orig, const std::string& relpath);
      explicit TrkCutHist(art::TFileDirectory tfdir);
    };

  } // namespace CutAndCount
} // namespace mu2e

#endif/*CutAndCountAnalysis_inc_TrkCutHist_hh*/
