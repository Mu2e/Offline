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
    {}
    StrawCluster::StrawCluster(ClusterType type,StrawId sid,
                               StrawEnd end,
                               double time,
                               float charge,
                               float ddist,
                               float phi,
                               float wdist,
                               float drifttime,
                               float proptime,
  			       art::Ptr<StrawGasStep> const& sgs,
			       art::Ptr<StepPointMC> const& stepmc,
			       CLHEP::HepLorentzVector const& cpos) : _type(type), _strawId(sid), _end(end), _time(time),
    _charge(charge), _ddist(ddist), _phi(phi),_wdist(wdist), _drifttime(drifttime), _proptime(proptime), _sgsptr(sgs), _spmcptr(stepmc), _cpos(cpos)
    {}

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
      << " drift distance " << _ddist
      << " phi " << _phi
      << " wire propagation distance " << _wdist 
      << *_sgsptr << std::endl;
    }
  }
}
