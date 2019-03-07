#ifndef TrackerConditions_Mu2eDetectorCache_hh
#define TrackerConditions_Mu2eDetectorCache_hh

//
// hold a set of run-dependent conditions objects
// and update them when needed
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/Mu2eDetectorMaker.hh"


namespace mu2e {
  class Mu2eDetectorCache : public ProditionsCache {
  public: 
    Mu2eDetectorCache(Mu2eDetectorConfig const& config):
      _name("Mu2eDetector"),_maker(config),
      _verbose(config.verbose()),_useDb(config.useDb()) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      if(!_mu2eMaterial_p) {
	_mu2eMaterial_p = std::make_unique<ProditionsHandle<Mu2eMaterial> >();
	_alignedDetector_p = std::make_unique<ProditionsHandle<Tracker> >();
      }

      auto& mu2eMaterial_h = *_mu2eMaterial_p;
      auto& alignedDetector_h = *_alignedDetector_p;

      // this pattern allows Mu2eDetector to depend on Mu2eMaterial
      // and AlignedTracker, assuming both have run dependence and
      // database dependence

      // get pointers to the currently valid objects
      auto mm = mu2eMaterial_h.getPtr(eid);
      auto tr = alignedDetector_h.getPtr(eid);
      
      // find overlapping valid ranges
      auto iov = mu2eMaterial_h.iov();
      iov.overlap(alignedDetector_h.iov());

      // this is the set of DB cid's they depend on
      ProditionsEntity::set_t s;
      s.merge(set_t(mm->getCids()));
      s.merge(set_t(tr->getCids()));
      // see if it is in cache
      auto p = find(s);
      if(!p) { // if not, make it
	if(_verbose>1) std::cout<< "making new Mu2eDetector " << std::endl;
	// since we do not have explicit db dependence
	// we can't differentiate on the useDb flag..
	p = _maker.fromFcl(mm,tr);
	p->addCids(s);
	push(p);
      } else {
	if(_verbose>1) std::cout<< "found Mu2eDetector in cache " << std::endl;
      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    Mu2eDetectorMaker _maker;
    int _verbose;
    bool _useDb;

    // this handle is not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<Mu2eMaterial> > _mu2eMaterial_p;
    std::unique_ptr<ProditionsHandle<Tracker> > _alignedDetector_p;

  };
};

#endif
