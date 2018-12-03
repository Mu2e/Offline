#ifndef DbExample_DetData3Cache_hh
#define DbExample_DetData3Cache_hh

#include <memory>
#include "DbExample/inc/DetData1.hh"
#include "DbExample/inc/DetData2.hh"
#include "DbExample/inc/ConditionsCache.hh"
#include "DbExample/inc/ConditionsHandle2.hh"
#include "DbExample/inc/DetData3Maker.hh"


namespace mu2e {
  class DetData3Cache : public ConditionsCache {
  public: 
    DetData3Cache():_name("DetData3") {}
    std::string const& name() const { return _name; }
    
    ConditionsCache::ret_t update(art::EventID const& eid) {
      
      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      // this has to be here because we can't create 
      // these handles in the constructor since this object
      // is created before the service is done creation.
      if(!_detData1_p) {
	_detData1_p = std::make_unique<ConditionsHandle2<DetData1> >();
	_detData2_p = std::make_unique<ConditionsHandle2<DetData2> >();
      }

      auto & _detData1_h = *_detData1_p;
      auto & _detData2_h = *_detData2_p;
      auto const& c1 = _detData1_h.get(eid);  // get the conditions objects
      auto const& c2 = _detData2_h.get(eid);  // this is dependant on
      auto iov = _detData1_h.iov();
      auto const& iov2 = _detData2_h.iov();
      iov.overlap(iov2);
      ConditionsEntity2::set_t s;
      s.insert(c1.getCids().begin(),c1.getCids().end());
      s.insert(c2.getCids().begin(),c2.getCids().end());
      auto p = find(s);
      if(!p) {
	std::cout<< "making new DetData3 " << std::endl;
	p = _maker.fromDb(c1,c2);
	p->addCids(s);
	push(p);
      } else {
	std::cout<< "found DetData3 " << std::endl;
      }
      return std::make_tuple(p,iov);
    }
    
  private:
    std::string _name;
    int _index;
    DetData3Maker _maker;
    std::unique_ptr<ConditionsHandle2<DetData1> > _detData1_p;
    std::unique_ptr<ConditionsHandle2<DetData2> > _detData2_p;
  };
};

#endif
