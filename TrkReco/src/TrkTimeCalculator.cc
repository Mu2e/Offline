//
// Object to calculate t0
//
// $Id: HelixFit.cc,v 1.12 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $
// $Date: 2014/07/10 14:47:26 $
//
// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
// root
#include "TH1F.h"
// C++
#include <vector>
#include <string>
#include <cmath>
using CLHEP::Hep3Vector;
using namespace std;
namespace mu2e
{

  TrkTimeCalculator::TrkTimeCalculator(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
//    _useflag(pset.get<std::vector<std::string>>("UseFlag")),
//    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier","Background"})),
    _avgDriftTime(pset.get<double>("AverageDriftTime",22.5)), 
    _useTOTdrift(pset.get<bool>("UseTOTDrift",true)),
    _beta(pset.get<double>("ParticleBeta",1.)),
    _shErr(pset.get<double>("StrawHitTimeErr",9.7)), // ns effective hit time res. without TOT
    _caloZOffset(pset.get<double>("CaloClusterZOffset",-120.0)), // WRT downstream face (mm)
    _caloT0Offset(pset.get<double>("TrkToCaloTimeOffset",-0.4)), // nanoseconds
    _caloT0Err(pset.get<double>("CaloTimeErr",0.5)) // nanoseconds
    { }

  TrkTimeCalculator::~TrkTimeCalculator() {}

  double TrkTimeCalculator::timeOfFlightTimeOffset(double hitz,double pitch) const {
    return hitz/(pitch*_beta*CLHEP::c_light);
  }

  double TrkTimeCalculator::comboHitTime(ComboHit const& ch,double pitch) {
    double tflt = timeOfFlightTimeOffset(ch.pos().z(),pitch);
    if (_useTOTdrift)
      return ch.correctedTime() - tflt; // use TOT to correct for drift
    else
      return ch.time() - tflt - _avgDriftTime; // otherwise make an average correction
  }

  double TrkTimeCalculator::caloClusterTime(CaloCluster const& cc,double pitch) const {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    Hep3Vector cog = ch->geomUtil().mu2eToTracker(ch->geomUtil().diskToMu2e( cc.diskId(), cc.cog3Vector())); 
    return cc.time() - timeOfFlightTimeOffset(cog.z()+_caloZOffset,pitch) + trkToCaloTimeOffset();
  }

}
