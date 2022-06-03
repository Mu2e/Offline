#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/DOCAStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXingUpdater.hh"
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
      config.bfcorr_ = fitconfig.bfieldCorr();
      config.ends_ = fitconfig.ends();
      config.tol_ = fitconfig.btol();
      // set the schedule for the meta-iterations
      std::vector<size_t> shualg;
      for(auto const& misetting : fitconfig.miConfig()) {
        MetaIterConfig mconfig(std::get<0>(misetting));
        config.schedule_.push_back(mconfig);
        shualg.push_back(std::get<1>(misetting));
      }
      // create the updaters requested
      unsigned ndoca(0); // indices into the updater config blocks
      unsigned nnull(0);
      auto const& nhusettings = fitconfig.nhuConfig();
      auto const& dhusettings = fitconfig.dhuConfig();
      auto const& sxusettings = fitconfig.sxuConfig();
      if(config.schedule_.size() != sxusettings.size())
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of KKStrawXing updaters" <<  std::endl;

      for( size_t imeta=0; imeta < config.schedule_.size(); ++imeta) {
        auto ialg = shualg[imeta];
        auto& miconfig = config.schedule_[imeta];
        if(ialg == StrawHitUpdaters::null) {
          auto const& nhusetting = nhusettings.at(nnull++);
          double maxdoca= std::get<0>(nhusetting);
          NullStrawHitUpdater nhupdater(maxdoca);
          miconfig.addUpdater(std::any(nhupdater));
        } else if(ialg == StrawHitUpdaters::DOCA) {
          auto const& dhusetting = dhusettings.at(ndoca++);
          double maxdoca= std::get<0>(dhusetting);
          double minddoca = std::get<1>(dhusetting);
          double maxddoca = std::get<2>(dhusetting);
          DOCAStrawHitUpdater dhupdater(maxdoca,minddoca,maxddoca);
          miconfig.addUpdater(std::any(dhupdater));
        } else {
          throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown updater " << ialg << std::endl;
        }
        //StrawXing updater too
        auto const& sxusetting = sxusettings.at(imeta);
        double maxdocasig= std::get<0>(sxusetting);
        double maxdoca = std::get<1>(sxusetting);
        double maxddoca = std::get<2>(sxusetting);
        KKStrawXingUpdater sxupdater(maxdocasig,maxdoca,maxddoca);
        miconfig.addUpdater(std::any(sxupdater));
      }
      // consistency test
      if(config.schedule_.size() != ndoca+nnull)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent StrawHitUpdater config "<< std::endl;
      return config;
    }
  }
}
