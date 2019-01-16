#ifndef TrackerConditions_StrawDriftCache_hh
#define TrackerConditions_StrawDriftCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/StrawDriftMaker.hh"


namespace mu2e {
  class StrawDriftCache : public ProditionsCache {
  public: 
    StrawDriftCache(StrawDriftConfig const& config):
      _name("StrawDrift"),_maker(config) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      //auto const& c1 = _c1_h.get(eid);
      DbIoV iov;
      iov.setMax();
      //iov = _c1_h.iov();
      ProditionsEntity::set_t s;
      //s.insert(_c1_h.cid());
      auto p = find(s);
      if(!p) {
	std::cout<< "making new StrawDrift " << std::endl;
	p = _maker.fromFcl();
	//p = _maker.fromDb(c1);
	//p->addCids(s);
	push(p);
      } else {
	std::cout<< "found StrawDrift in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    StrawDriftMaker _maker;

  };
};

#endif
