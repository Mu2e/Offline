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
      std::vector<size_t> shualg;
      for(auto const& misetting : fitconfig.miConfig()) {
        MetaIterConfig mconfig(std::get<0>(misetting));
        config.schedule_.push_back(mconfig);
        shualg.push_back(std::get<1>(misetting));
      }
      // create the updaters requested
      unsigned ndoca(0);
      auto const& dhusettings = fitconfig.dhuConfig();
      for( size_t imeta=0; imeta < config.schedule_.size(); ++imeta) {
        auto ialg = shualg[imeta];
        auto& miconfig = config.schedule_[imeta];
        if(ialg == StrawHitUpdater::null) {
          miconfig.addUpdater(std::any(NullStrawHitUpdater()));
        } else if(ialg == StrawHitUpdater::DOCA) {
          auto const& dhusetting = dhusettings.at(ndoca++);
          double maxdoca= std::get<0>(dhusetting);
          double minddoca = std::get<1>(dhusetting);
          double maxddoca = std::get<2>(dhusetting);
          DOCAStrawHitUpdater dhupdater(maxdoca,minddoca,maxddoca);
          miconfig.addUpdater(std::any(dhupdater));
        } else {
          throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown updater " << ialg << std::endl;
        }
      }
      return config;
    }
  }
}
