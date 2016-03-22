// Andrei Gaponenko, 2016

#include "CutAndCountAnalysis/inc/TrkCutHist.hh"

namespace mu2e {
  namespace CutAndCount {

    art::TFileDirectory TrkCutHist::getdir(art::TFileDirectory orig, const std::string& relpath) {
      return relpath.empty() ? orig : orig.mkdir(relpath);
    }

    TrkCutHist::TrkCutHist(art::TFileDirectory tfdir, const std::string& relpath)
      : TrkCutHist(getdir(tfdir, relpath))
    {}

    TrkCutHist::TrkCutHist(art::TFileDirectory tf)
      : trkqual{tf.make<TH1D>("trkqual", "trkqual before cut", 100, 0., 1.)}
      , td{tf.make<TH1D>("td", "Track tan(lambda) before cut", 100, 0.5, 1.5)}
      , d0{tf.make<TH1D>("d0", "Track d0 before cut", 300, -150., +150.)}
      , rmax{tf.make<TH1D>("rmax", "Track d0+2/om  before cut", 120, 300., 900.)}
      , t0{tf.make<TH1D>("t0", "Track t0  before cut", 170, 0., 1700.)}
      , caloMatchChi2{tf.make<TH1D>("caloMatchCHi2", "Calo match chi2 before cut", 100, 0., 300.)}
      , caloClusterEnergy{tf.make<TH1D>("caloClusterEnergy", "Calo cluster energy before cut", 150, 0., 150.)}
      , momentum{tf.make<TH1D>("momentum", "Track momentum  before cut", 500, 98., 108.)}
    {
      trkqual->Sumw2();
      td->Sumw2();
      d0->Sumw2();
      rmax->Sumw2();
      t0->Sumw2();
      caloMatchChi2->Sumw2();
      caloClusterEnergy->Sumw2();
      momentum->Sumw2();
    }
  } // namespace CutAndCount
} // namespace mu2e
