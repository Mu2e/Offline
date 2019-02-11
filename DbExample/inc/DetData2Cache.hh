#ifndef DbExample_DetData2Cache_hh
#define DbExample_DetData2Cache_hh

#include "DbExample/inc/ConditionsCache.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbExample/inc/DetData2Maker.hh"


namespace mu2e {
  class DetData2Cache : public ConditionsCache {
  public: 
    DetData2Cache():_name("DetData2") {}
    std::string const& name() const { return _name; }

    ConditionsCache::ret_t update(art::EventID const& eid) {
      
      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      auto const& c1 = _c1_h.get(eid);
      auto const& c2 = _c2_h.get(eid);
      auto iov = _c1_h.iov();
      auto const& iov2 = _c2_h.iov();
      iov.overlap(iov2);
      ConditionsEntity2::set_t s;
      s.insert(_c1_h.cid());
      s.insert(_c2_h.cid());
      auto p = find(s);
      if(!p) {
	std::cout<< "making new DetData2 " << std::endl;
	p = _maker.fromDb(c1,c2);
	p->addCids(s);
	push(p);
      } else {
	std::cout<< "found DetData2 " << std::endl;
      }
      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    DetData2Maker _maker;
    DbHandle<TstCalib1> _c1_h;
    DbHandle<TstCalib2> _c2_h;
  };
};

#endif
