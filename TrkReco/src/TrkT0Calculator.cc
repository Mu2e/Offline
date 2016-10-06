//
// Object to calculate t0
//
// $Id: HelixFit.cc,v 1.12 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkT0Calculator.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
// root
#include "TH1F.h"
// C++
#include <vector>
#include <string>
#include <math.h>
#include <cmath>
using CLHEP::Hep3Vector;
using namespace std;
using namespace boost::accumulators;
namespace mu2e 
{

  TrkT0Calculator::TrkT0Calculator(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _useflag(pset.get<std::vector<std::string>>("UseFlag")),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier","DeltaRay","Isolated"}))
    {} 

  TrkT0Calculator::~TrkT0Calculator() {}

  void TrkT0Calculator::updateT0(TimeCluster& tc, StrawHitCollection const& shcol, TrkFitDirection const& fdir){


  }
  void TrkT0Calculator::updateT0(HelixSeed& hs, StrawHitCollection const& shcol, TrkFitDirection const& fdir) {

  }

}
 
