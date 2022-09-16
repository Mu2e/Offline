#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
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
      config.pdchisq_ = fitconfig.dparams();
      config.divgap_ = fitconfig.dgap();
      config.bfcorr_ = fitconfig.bfieldCorr();
      config.ends_ = fitconfig.ends();
      config.tol_ = fitconfig.btol();
      // set the schedule for the meta-iterations
      std::vector<int> shualg;
      for(auto const& misetting : fitconfig.miConfig()) {
        MetaIterConfig mconfig(std::get<0>(misetting));
        config.schedule_.push_back(mconfig);
        auto alg = StrawHitUpdaters::algo(std::get<1>(misetting));
        if(alg == StrawHitUpdaters::unknown)
          throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown StrawHitUpdater " << std::get<1>(misetting) << std::endl;
        shualg.push_back(static_cast<int>(alg));
      }
      // create the updaters requested
      unsigned nptca, nnull, nann, nbkg, ncomb, nnone;
      nptca = nnull = nann = nbkg = ncomb = nnone= 0; // count how many updater configs have been seen
      std::vector<CAStrawHitUpdater::CASHUConfig> cashusettings;
      std::vector<ANNStrawHitUpdater::ANNSHUConfig> annshusettings;
      std::vector<BkgStrawHitUpdater::BkgSHUConfig> bkgshusettings;
      std::vector<CombinatoricStrawHitUpdater::CSHUConfig> chusettings;
      // specific updaters can be empty, so fetch config data with a default empty vector
      cashusettings = fitconfig.cashuConfig().value_or(cashusettings);
      annshusettings = fitconfig.annshuConfig().value_or(annshusettings);
      bkgshusettings = fitconfig.bkgshuConfig().value_or(bkgshusettings);
      chusettings = fitconfig.chuConfig().value_or(chusettings);
      auto const& sxusettings = fitconfig.sxuConfig();
      if(config.schedule_.size() != sxusettings.size())
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of KKStrawXing updaters" <<  std::endl;
      for( size_t imeta=0; imeta < config.schedule_.size(); ++imeta) {
        auto alg = shualg[imeta];
        auto& miconfig = config.schedule_[imeta];
        if(alg == StrawHitUpdaters::CA) {
          miconfig.addUpdater(std::any(CAStrawHitUpdater(cashusettings.at(nptca++))));
        } else if(alg == StrawHitUpdaters::ANN) {
          miconfig.addUpdater(std::any(ANNStrawHitUpdater(annshusettings.at(nann++))));
        } else if(alg == StrawHitUpdaters::Bkg) {
          miconfig.addUpdater(std::any(BkgStrawHitUpdater(bkgshusettings.at(nbkg++))));
        } else if(alg == StrawHitUpdaters::Combinatoric) {
          miconfig.addUpdater(std::any(CombinatoricStrawHitUpdater(chusettings.at(ncomb++))));
        } else if(alg == StrawHitUpdaters::none) {
          ++nnone;
        } else {
          throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown StrawHitUpdater " << alg << std::endl;
        }
        //StrawXing updater too; these must always be present
        miconfig.addUpdater(std::any(StrawXingUpdater(sxusettings.at(imeta))));
      }
      // consistency test
      if(config.schedule_.size() != nptca+nnull+nann+nbkg+ncomb+nnone)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent StrawHitUpdater config "<< std::endl;
      return config;
    }
  }
}
