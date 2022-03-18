#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitUpdater.hh"
#include "KinKal/Detector/WireHitStructs.hh"
#include <iostream>

namespace mu2e {
  using KinKal::Config;
  using KinKal::MetaIterConfig;
  namespace Mu2eKinKal {
    Config makeConfig(KinKalConfig const& fitconfig) {
      Config config;
      // fill configuration
      config.plevel_ = static_cast<Config::printLevel>(fitconfig.printLevel());
      config.minndof_ = fitconfig.minndof();
      config.maxniter_ = fitconfig.maxniter();
      config.dwt_ = fitconfig.dwt();
      config.convdchisq_ = fitconfig.convdchisq();
      config.divdchisq_ = fitconfig.divdchisq();
      config.pdchi2_ = fitconfig.dparams();
      config.tbuff_ = fitconfig.tBuffer();
      config.bfcorr_ = fitconfig.bfieldCorr();
      config.tol_ = fitconfig.btol();
      // set the schedule for the meta-iterations
      for(auto const& misetting : fitconfig.miConfig()) {
        MetaIterConfig mconfig(std::get<0>(misetting));
        config.schedule_.push_back(mconfig);
      }
      // create the updaters requested; these don't have to exist, but if they
      // do, their dimension must match the meta-iterations
      if(fitconfig.shuConfig().size() >0 && fitconfig.shuConfig().size() != config.schedule_.size())
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: KKStrawHitUpdaters don't match meta-iterations" << std::endl;
      size_t imeta(0);
      for(auto const& shusetting : fitconfig.shuConfig()) {
        double maxdoca= std::get<0>(shusetting);
        double minprob = std::get<1>(shusetting);
        double minddoca = std::get<2>(shusetting);
        double maxddoca = std::get<3>(shusetting);
        KKStrawHitUpdater shupdater(maxdoca,minprob,minddoca,maxddoca);
        config.schedule_[imeta++].addUpdater(std::any(shupdater));
      }
      return config;
    }
  }
}
