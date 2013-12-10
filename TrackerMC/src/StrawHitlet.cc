//
// StrawHitlet
// $Id: StrawHitlet.cc,v 1.2 2013/12/10 01:32:51 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/10 01:32:51 $
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawHitlet.hh"
// general includes

using namespace std;
namespace mu2e {

  StrawHitlet::StrawHitlet() : _type(unknown), _strawIndex(0), _end(StrawEnd::unknown), _time(0.0), _charge(0.0), _wdist(0.0)
  {}

  StrawHitlet::StrawHitlet(HitletType type,StrawIndex sindex,
      StrawEnd end,
      double time,
      double charge,
      double wdist,
      art::Ptr<StepPointMC> const& stepmc) : _type(type), _strawIndex(sindex), _end(end), _time(time),
  _charge(charge), _wdist(wdist), _stepMC(stepmc) 
  {}

  StrawHitlet::StrawHitlet(const StrawHitlet& other) :
    _type(other._type), _strawIndex(other._strawIndex), _end(other._end),
    _time(other._time), _charge(other._charge), _wdist(other._wdist) 
  {}

// delegating constructors in C++11!
  StrawHitlet::StrawHitlet(const StrawHitlet& primary, StrawIndex const& index, double xfactor) :
  StrawHitlet(primary) {
    _type = xtalk;
    _strawIndex = index;
    _charge *= xfactor;
  }

  StrawHitlet::StrawHitlet(const StrawHitlet& primary, double deltat) : StrawHitlet(primary) {
    _time += deltat;
  }

  StrawHitlet& StrawHitlet::operator=(StrawHitlet const& other ) {
    if(&other != this){
      _type = other._type;
      _strawIndex = other._strawIndex;
      _end = other._end;
      _time = other._time;
      _charge = other._charge;
      _wdist = other._wdist;
      _stepMC = other._stepMC;
    }
    return *this;
  }

  void StrawHitlet::print(std::ostream& ost, bool doEndl) const {
    ost << "StrawHitlet of type " << _type
    << " for straw index " << _strawIndex
    << " end " << _end 
    << " time " << _time
    << " charge " << _charge
    << " wire distance " << _wdist
    << " StepPointMC ";
    _stepMC->print(ost,doEndl);
  }

}

