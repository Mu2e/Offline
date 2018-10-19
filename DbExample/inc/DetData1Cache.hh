#ifndef DbExample_DetData1Cache_hh
#define DbExample_DetData1Cache_hh

#include "DbExample/inc/ConditionsCache.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbExample/inc/DetData1Maker.hh"


namespace mu2e {
  class DetData1Cache : public ConditionsCache {
  public: 
    DetData1Cache():_name("DetData1") {}
    std::string const& name() const { return _name; }

    ConditionsCache::ret_t update(art::EventID const& eid) {

      auto const& c1 = _c1_h.get(eid);
      auto iov = _c1_h.iov();
      ConditionsCache::set_t s;
      s.insert(_c1_h.cid());
      auto p = find(s);
      if(!p) {
	p = _maker.fromDb(c1);
	std::cout<< "making new DetData1 " << std::endl;
      } else {
	std::cout<< "found DetData1 " << std::endl;
      }
      push(p,s);
      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    DetData1Maker _maker;
    DbHandle<TstCalib1> _c1_h;
  };
};

#endif
