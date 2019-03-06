#ifndef TrackerConditions_StrawPhysicsCache_hh
#define TrackerConditions_StrawPhysicsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawPhysicsMaker.hh"


namespace mu2e {
  class StrawPhysicsCache : public ProditionsCache {
  public: 
    StrawPhysicsCache(StrawPhysicsConfig const& config):
      _name("StrawPhysics"),_maker(config),
      _verbose(config.verbose()),_useDb(config.useDb()) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      // these have to be created on first use because we can't 
      // create them in the constructor since this object
      // is created before the service is done creation.
       if(!_strawDrift_p) {
	_strawDrift_p = std::make_unique<ProditionsHandle<StrawDrift> >();
      }
      auto & strawDrift_h = *_strawDrift_p;

      auto sd = strawDrift_h.getPtr(eid);
      auto iov = strawDrift_h.iov();
      ProditionsEntity::set_t s;
      s.merge(set_t(sd->getCids()));
      auto p = find(s);
      if(!p) {
	if(_verbose>1) std::cout<< "making new StrawPhysics " << std::endl;
	p = _maker.fromFcl(sd);
	//p = _maker.fromDb(sd);
	p->addCids(s);
	push(p);
      } else {
	if(_verbose>1) std::cout<< "found StrawPhysics in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    StrawPhysicsMaker _maker;
    int _verbose;
    bool _useDb;

    // this handle is not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<StrawDrift> > _strawDrift_p;
  };
};

#endif
