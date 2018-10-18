//
// StrawCluster
//
// Original author David Brown, LBNL
//
// mu2e includes
#include "TrackerMC/inc/StrawCluster.hh"
// general includes

using namespace std;
namespace mu2e {
  namespace TrackerMC {
    StrawCluster::StrawCluster() : _type(unknown), _strawId(0), _end(StrawEnd::cal), _time(0.0), _charge(0.0), _ddist(0.0),_phi(0.0),  _wdist(0.0), _drifttime(0.0), _proptime(0.0)
    {}//JB: added phi
    StrawCluster::StrawCluster(ClusterType type,StrawId sid,
                               StrawEnd end,
                               double time,
                               double charge,
                               double ddist,
                               double phi,//JB: added
                               double wdist,
                               double drifttime,
                               double proptime,
                               art::Ptr<StepPointMC> const& stepmc,
                               CLHEP::HepLorentzVector const& cpos) : _type(type), _strawId(sid), _end(end), _time(time),
    _charge(charge), _ddist(ddist), _phi(phi),_wdist(wdist), _drifttime(drifttime), _proptime(proptime), _stepMC(stepmc), _cpos(cpos)
    {}//JB: added phi

    StrawCluster::StrawCluster(const StrawCluster& other) :
    _type(other._type), _strawId(other._strawId), _end(other._end),
    _time(other._time), _charge(other._charge), _ddist(other._ddist), _phi(other._phi), _wdist(other._wdist), _drifttime(other._drifttime), _proptime(other._proptime), _stepMC(other._stepMC), _cpos(other._cpos)
    {}

    // delegating constructors in C++11!
    StrawCluster::StrawCluster(const StrawCluster& primary, StrawId const& id, double xfactor) :
    StrawCluster(primary) {
      _type = xtalk;
      _strawId = id;
      _charge *= xfactor;
    }

    StrawCluster::StrawCluster(const StrawCluster& primary, double deltat) : StrawCluster(primary) {
      _time += deltat;
    }

    StrawCluster& StrawCluster::operator=(StrawCluster const& other ) {
      if(&other != this){
        _type = other._type;
        _strawId = other._strawId;
        _end = other._end;
        _time = other._time;
        _charge = other._charge;
        _ddist = other._ddist;
        _phi = other._phi; //JB: added phi
        _wdist = other._wdist;
        _drifttime = other._drifttime;
        _proptime = other._proptime;
        _stepMC = other._stepMC;
        _cpos = other._cpos;
      }
      return *this;
    }

    void StrawCluster::print(std::ostream& ost, bool doEndl) const {
      ost << "StrawCluster of type " << _type
      << " for straw id " << _strawId
      << " end " << _end
      << " time " << _time
      << " charge " << _charge
      << " drift distance " << _ddist
      << " phi " << _phi
      << " wire propagation distance " << _wdist
      << " StepPointMC ";
      _stepMC->print(ost,doEndl); //JB: added phi
    }
  }
}
