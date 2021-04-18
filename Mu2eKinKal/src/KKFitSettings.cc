#include "Mu2eKinKal/inc/KKFitSettings.hh"
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
      for(auto const& misetting : fitconfig.mconfig()) {
	MetaIterConfig mconfig;
	mconfig.temp_ = std::get<0>(misetting);
	mconfig.convdchisq_ = std::get<1>(misetting);
	mconfig.divdchisq_ = std::get<2>(misetting);
// setup the updaters as appropriate
	int updateflag = std::get<3>(misetting);
	if(updateflag != 0){
	  std::cout << "updateflag = " << updateflag << std::endl;
	}
	mconfig.miter_ = nmiter++;
	config.schedule_.push_back(mconfig);
      }
    // set the hit updating
//      KKStrawHitUpdater shupdater(std::get<0>(shusetting), std::get<1>(shusetting), kkfit_.nullDimension());
//      unsigned minmeta = std::get<2>(shusetting);
//      unsigned maxmeta = std::get<3>(shusetting);
//      unsigned maxmeta = std::get<4>(shusetting);
//      if(maxmeta < minmeta || schedule.size() < maxmeta)
//	throw cet::exception("RECO")<<"mu2e::LoopHelixFit: Hit updater configuration error"<< endl;
//      for(unsigned imeta=minmeta; imeta<=maxmeta; imeta++)
//	schedule[imeta].updaters_.push_back(shupdater);
//    }


      return config;
    }
  }
}
