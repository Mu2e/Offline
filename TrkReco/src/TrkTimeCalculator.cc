#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

namespace mu2e
{

  TrkTimeCalculator::TrkTimeCalculator(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    // _useflag(pset.get<std::vector<std::string>>("UseFlag")),
    // _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier","Background"})),
    _avgDriftTime(pset.get<double>("AverageDriftTime",22.5)),
    _useTOTdrift(pset.get<bool>("UseTOTDrift",true)),
    _beta(pset.get<double>("StrawHitBeta",1.)),
    _shErr(pset.get<double>("StrawHitTimeErr",9.7)), // ns effective hit time res. without TOT
    _caloZOffset(pset.get<double>("CaloClusterZOffset",-120.0)), // WRT downstream face (mm)
    _caloTimeOffset(pset.get<double>("TrkToCaloTimeOffset",-0.4)), // nanoseconds
    _caloTimeErr(pset.get<double>("CaloTimeErr",0.5)) // nanoseconds
  { }

  TrkTimeCalculator::TrkTimeCalculator(const Config& config) :
    _debug(         config.debug()),
    //_useflag(       config.useFlag()),
    //_dontuseflag(   config.dontUseFlag()),
    _avgDriftTime(  config.avgDriftTime()),
    _useTOTdrift(   config.useTOTdrift()),
    _beta(          config.beta()),
    _shErr(         config.shErr()),
    _caloZOffset(   config.caloZOffset()),
    _caloTimeOffset(  config.caloTimeOffset()),
    _caloTimeErr(     config.caloTimeErr())
  {}


  double TrkTimeCalculator::timeOfFlightTimeOffset(double hitz,double pitch) const
  {
    return hitz/(pitch*_beta*CLHEP::c_light);
  }

  double TrkTimeCalculator::comboHitTime(ComboHit const& ch,double pitch)
  {
    double tflt = timeOfFlightTimeOffset(ch.pos().z(),pitch);
    if (_useTOTdrift)
      return ch.correctedTime() - tflt; // use TOT to correct for drift
    else
      return ch.time() - tflt - _avgDriftTime; // otherwise make an average correction
  }

  double TrkTimeCalculator::caloClusterTime(CaloCluster const& cc,double pitch) const
  {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    CLHEP::Hep3Vector cog = ch->geomUtil().mu2eToTracker(ch->geomUtil().diskToMu2e( cc.diskID(), cc.cog3Vector()));
    return cc.time() - timeOfFlightTimeOffset(cog.z()+_caloZOffset,pitch) + trkToCaloTimeOffset();
  }

}
