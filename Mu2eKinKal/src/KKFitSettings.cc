#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/ANNStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/BkgStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/StrawXingUpdater.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
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
        shualg.push_back(std::get<1>(misetting));
      }
      // create the updaters requested
      unsigned nptca, nnull, nann, nbkg, ncomb, nnone;
      nptca = nnull = nann = nbkg = ncomb = nnone= 0; // count how many updater configs have been seen
      std::vector<NullStrawHitUpdater::NSHUConfig> nhusettings;
      std::vector<PTCAStrawHitUpdater::PSHUConfig> pshusettings;
      std::vector<ANNStrawHitUpdater::ANNSHUConfig> annshusettings;
      std::vector<BkgStrawHitUpdater::BkgSHUConfig> bkgshusettings;
      std::vector<CombinatoricStrawHitUpdater::CSHUConfig> chusettings;
      // specific updaters can be empty, so fetch config data with a default empty vector
      nhusettings = fitconfig.nhuConfig().value_or(nhusettings);
      pshusettings = fitconfig.pshuConfig().value_or(pshusettings);
      annshusettings = fitconfig.annshuConfig().value_or(annshusettings);
      bkgshusettings = fitconfig.bkgshuConfig().value_or(bkgshusettings);
      chusettings = fitconfig.chuConfig().value_or(chusettings);
      auto const& sxusettings = fitconfig.sxuConfig();
      if(config.schedule_.size() != sxusettings.size())
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of KKStrawXing updaters" <<  std::endl;
      for( size_t imeta=0; imeta < config.schedule_.size(); ++imeta) {
        auto ialg = shualg[imeta];
        auto& miconfig = config.schedule_[imeta];
        if(ialg == StrawHitUpdaters::null) {
          miconfig.addUpdater(std::any(NullStrawHitUpdater(nhusettings.at(nnull++))));
        } else if(ialg == StrawHitUpdaters::PTCA) {
          miconfig.addUpdater(std::any(PTCAStrawHitUpdater(pshusettings.at(nptca++))));
        } else if(ialg == StrawHitUpdaters::ANN) {
          miconfig.addUpdater(std::any(ANNStrawHitUpdater(annshusettings.at(nann++))));
        } else if(ialg == StrawHitUpdaters::Bkg) {
          miconfig.addUpdater(std::any(BkgStrawHitUpdater(bkgshusettings.at(nbkg++))));
        } else if(ialg == StrawHitUpdaters::Combinatoric) {
          miconfig.addUpdater(std::any(CombinatoricStrawHitUpdater(chusettings.at(ncomb++))));
        } else if(ialg == StrawHitUpdaters::none) {
          ++nnone;
        } else {
          throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown updater " << ialg << std::endl;
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
