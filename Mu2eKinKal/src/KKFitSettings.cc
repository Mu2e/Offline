#include "Mu2eKinKal/inc/KKFitSettings.hh"
#include "Mu2eKinKal/inc/KKStrawHitUpdater.hh"
#include "KinKal/Detector/WireHitStructs.hh"
#include <iostream>

namespace mu2e {
  using KinKal::Config;
  using KinKal::MetaIterConfig;
  namespace Mu2eKinKal {
    Config makeConfig(KKConfig const& fitconfig) {
      Config config;
      // fill configuration
      config.maxniter_ = fitconfig.maxniter();
      config.dwt_ = fitconfig.dwt();
      config.pdchi2_ = fitconfig.dparams();
      config.tbuff_ = fitconfig.tBuffer();
      config.tol_ = fitconfig.btol();
      config.minndof_ = fitconfig.minndof();
      config.bfcorr_ = static_cast<Config::BFCorr>(fitconfig.bfieldCorr());
      config.plevel_ = static_cast<Config::printLevel>(fitconfig.printLevel());
      // set the schedule for the meta-iterations
      unsigned nmiter(0);
      for(auto const& misetting : fitconfig.miConfig()) {
	MetaIterConfig mconfig;
	mconfig.miter_ = nmiter++;
	mconfig.temp_ = std::get<0>(misetting);
	mconfig.convdchisq_ = std::get<1>(misetting);
	mconfig.divdchisq_ = std::get<2>(misetting);
	config.schedule_.push_back(mconfig);
      }
// create the updaters requested
      for(auto const& shusetting : fitconfig.shuConfig()) {
	size_t imeta = std::get<4>(shusetting);
	if(imeta >= config.schedule_.size())
	  throw cet::exception("RECO")<<"mu2e::KKFitSettings: Invalid Meta-iteration for KKStrawHitUpdater" << std::endl;
	double mindoca = std::get<0>(shusetting);
	double maxdoca= std::get<1>(shusetting);
	double maxchi = std::get<2>(shusetting);
	WireHitState::Dimension ndim = static_cast<WireHitState::Dimension>(std::get<3>(shusetting));
	KKStrawHitUpdater shupdater(mindoca,maxdoca,maxchi,ndim);
	config.schedule_[imeta].updaters_.push_back(shupdater);
      }

      return config;
    }
  }
}
