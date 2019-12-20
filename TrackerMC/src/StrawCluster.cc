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
    StrawCluster::StrawCluster() : _type(unknown), _strawId(0), _end(StrawEnd::cal), _time(0.0), _charge(0.0),  _drifttime(0.0), _proptime(0.0)
    {}

    StrawCluster::StrawCluster(ClusterType type,StrawId sid,
                               StrawEnd end,
                               float time,
                               float charge,
			       float wdist,
			       StrawPosition const& pos,
                               float drifttime,
                               float proptime,
			       art::Ptr<StrawGasStep> const& sgsptr,
			       float ctime) : _type(type), _strawId(sid), _end(end), _time(time),
    _charge(charge), _wdist(wdist), _pos(pos), _drifttime(drifttime), _proptime(proptime), _sgsptr(sgsptr), _ctime(ctime) {}

    // delegating constructors in C++11!
    StrawCluster::StrawCluster(const StrawCluster& primary, StrawId const& id, float xfactor) :
    StrawCluster(primary) {
      _type = xtalk;
      _strawId = id;
      _charge *= xfactor;
    }

    StrawCluster::StrawCluster(const StrawCluster& primary, double deltat) : StrawCluster(primary) {
      _time += deltat;
    }

    void StrawCluster::print(std::ostream& ost, bool doEndl) const {
      ost << "StrawCluster of type " << _type
      << " for straw id " << _strawId
      << " end " << _end
      << " time " << _time
      << " charge " << _charge
      << " drift distance " << _pos.Rho()
      << " phi " << _pos.Phi()
      << " wire propagation distance " << _wdist
      << " " << *_sgsptr << std::endl;
    }
  }
}
