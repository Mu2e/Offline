#ifndef TrackerConditions_StrawElectronicsCache_hh
#define TrackerConditions_StrawElectronicsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "TrackerConditions/inc/StrawElectronicsMaker.hh"


namespace mu2e {
  class StrawElectronicsCache : public ProditionsCache {
  public: 
    StrawElectronicsCache(StrawElectronicsConfig const& config):
      _name("StrawElectronics"),_maker(config) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      ProditionsEntity::set_t s;
      DbIoV iov;
      iov.setMax();
      auto p = find(s);
      if(!p) {
	std::cout<< "making new StrawElectronics " << std::endl;
	p = _maker.fromFcl();
	p->addCids(s);
	push(p);
      } else {
	std::cout<< "found StrawElectronics in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    StrawElectronicsMaker _maker;

  };
};

#endif
