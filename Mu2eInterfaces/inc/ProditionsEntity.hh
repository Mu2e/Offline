#ifndef Mu2eInterfaces_ProditionsEntity_hh
#define Mu2eInterfaces_ProditionsEntity_hh

#include <set>
#include <memory>

namespace mu2e {
  class ProditionsEntity {
  public:
    typedef std::shared_ptr<ProditionsEntity> ptr;
    typedef std::set<int> set_t;

    virtual ~ProditionsEntity() {}
    virtual std::string const& name() const =0 ;
    set_t const& getCids() const { return _cids; }
    void addCids(set_t const& s) { _cids.insert(s.begin(),s.end()); }
    virtual void print( std::ostream& ) const {}
  private:
    set_t _cids;
  };
};

#endif
