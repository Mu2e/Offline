#ifndef TrackerConditions_Mu2eMaterialCache_hh
#define TrackerConditions_Mu2eMaterialCache_hh

//
// This Proditions entitiy cache shoudl only ever hold one
// because Mu2eMaterail has pointers to BTrk singletons
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/Mu2eMaterialMaker.hh"


namespace mu2e {
  class Mu2eMaterialCache : public ProditionsCache {
  public: 
    Mu2eMaterialCache(Mu2eMaterialConfig const& config):
      _name("Mu2eMaterial"),_maker(config),
      _verbose(config.verbose()) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      // this pattern only allows one copy
      DbIoV iov;
      iov.setMax(); // all runs
      ProditionsEntity::set_t s;
      auto p = find(s);
      if(!p) {
	if(_verbose>1) std::cout<< "making new Mu2eMaterial " << std::endl;
	p = _maker.fromFcl();
	p->addCids(s);
	push(p);

	if(_verbose>2) p->print(std::cout);

      } else {
	if(_verbose>1) std::cout<< "found Mu2eMaterial in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    Mu2eMaterialMaker _maker;
    int _verbose;

  };
};

#endif
