#include "Offline/Mu2eKinKal/inc/KKFitSettings.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/GeneralUtilities/inc/splitLine.hh"
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
      // create the updaters requested
      std::vector<CADSHU::Config> cadshusettings;
      std::vector<DriftANNSHU::Config> driftannshusettings;
      std::vector<BkgANNSHU::Config> bkgannshusettings;
      std::vector<Chi2SHU::Config> chi2shusettings;
      // specific updaters can be empty, so fetch config data with a default empty vector
      cadshusettings = fitconfig.cadshuConfig().value_or(cadshusettings);
      driftannshusettings = fitconfig.annshuConfig().value_or(driftannshusettings);
      bkgannshusettings = fitconfig.bkgshuConfig().value_or(bkgannshusettings);
      chi2shusettings = fitconfig.combishuConfig().value_or(chi2shusettings);
      // straw material updater must always be here
      auto const& sxusettings = fitconfig.sxuConfig();
      // set the schedule for the meta-iterations
      unsigned ncadshu(0), nann(0), nbkg(0), ncomb(0), nnone(0), nsxu(0);
      for(auto const& misetting : fitconfig.miConfig()) {
        MetaIterConfig miconfig(std::get<0>(misetting));
        // parse StrawHit updaters, and add to the config of this meta-iteraion
        std::vector<std::string> anames;
        splitLine( std::get<1>(misetting), ":", anames);
        for(auto const& aname : anames) {
          auto alg = StrawHitUpdaters::algo(aname);
          if(alg == StrawHitUpdaters::CAD) {
            miconfig.addUpdater(std::any(CADSHU(cadshusettings.at(ncadshu++))));
          } else if(alg == StrawHitUpdaters::DriftANN) {
            miconfig.addUpdater(std::any(DriftANNSHU(driftannshusettings.at(nann++))));
          } else if(alg == StrawHitUpdaters::BkgANN) {
            miconfig.addUpdater(std::any(BkgANNSHU(bkgannshusettings.at(nbkg++))));
          } else if(alg == StrawHitUpdaters::Chi2) {
            miconfig.addUpdater(std::any(Chi2SHU(chi2shusettings.at(ncomb++))));
          } else if(alg == StrawHitUpdaters::none) {
            ++nnone;
          } else {
            throw cet::exception("RECO")<<"mu2e::KKFitSettings: unknown StrawHitUpdater " << alg << std::endl;
          }
        }
        // pad straw xing updaters if necessary
        miconfig.addUpdater(std::any(StrawXingUpdater(sxusettings.at(nsxu))));
        if(sxusettings.size()> nsxu+1)nsxu++;
        config.schedule_.push_back(miconfig);
      }
      // consistency checks
      if(cadshusettings.size() != ncadshu)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of CA StrawHit updaters" <<  std::endl;
      if(driftannshusettings.size() != nann)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of ANN StrawHit updaters" <<  std::endl;
      if(bkgannshusettings.size() != nbkg)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of Bkg StrawHit updaters" <<  std::endl;
      if(chi2shusettings.size() != ncomb)
        throw cet::exception("RECO")<<"mu2e::KKFitSettings: inconsistent number of Combi StrawHit updaters" <<  std::endl;
      return config;
    }
  }
}
