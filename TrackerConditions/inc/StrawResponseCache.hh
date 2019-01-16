#ifndef TrackerConditions_StrawResponseCache_hh
#define TrackerConditions_StrawResponseCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/StrawResponseMaker.hh"


namespace mu2e {
  class StrawResponseCache : public ProditionsCache {
  public: 
    StrawResponseCache(StrawResponseConfig const& config):
      _name("StrawResponse"),_maker(config) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      if(!_strawDrift_p) {
	_strawDrift_p = std::make_unique<ProditionsHandle<StrawDrift> >();
	_strawElectronics_p = std::make_unique<ProditionsHandle<StrawElectronics> >();
	_strawPhysics_p = std::make_unique<ProditionsHandle<StrawPhysics> >();
      }
      auto & strawDrift_h = *_strawDrift_p;
      auto & strawElectronics_h = *_strawElectronics_p;
      auto & strawPhysics_h = *_strawPhysics_p;

      // find the smallest common IOV range
      auto sd = strawDrift_h.getPtr(eid);
      auto iov = strawDrift_h.iov();
      ProditionsEntity::set_t s = sd->getCids();

      auto se = strawElectronics_h.getPtr(eid);
      iov.overlap(strawDrift_h.iov());
      s.merge(set_t(se->getCids()));

      auto sp = strawPhysics_h.getPtr(eid);
      iov.overlap(strawPhysics_h.iov());
      s.merge(set_t(sp->getCids()));

      auto p = find(s);
      if(!p) {
	std::cout<< "making new StrawResponse " << std::endl;
	p = _maker.fromFcl(sd,se,sp);
	//p = _maker.fromDb(..);
	p->addCids(s);
	push(p);
      } else {
	std::cout<< "found StrawResponse in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    StrawResponseMaker _maker;
    std::unique_ptr<ProditionsHandle<StrawDrift> > _strawDrift_p;
    std::unique_ptr<ProditionsHandle<StrawElectronics> > _strawElectronics_p;
    std::unique_ptr<ProditionsHandle<StrawPhysics> > _strawPhysics_p;
  };
};

#endif
