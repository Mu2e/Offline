#ifndef Mu2eInterfaces_ProditionsEntity_hh
#define Mu2eInterfaces_ProditionsEntity_hh

#include <set>
#include <memory>

namespace mu2e {
  class ProditionsEntity {
  public:
    typedef std::shared_ptr<ProditionsEntity> ptr;
    typedef std::set<int> set_t;

    ProditionsEntity(std::string const& name):_name(name) {}
    virtual ~ProditionsEntity() = default;
    ProditionsEntity( ProditionsEntity const&  ) = default;
    ProditionsEntity( ProditionsEntity&&       ) = default;
    ProditionsEntity& operator=(ProditionsEntity const&  ) = delete;
    ProditionsEntity& operator=(ProditionsEntity &&      ) = delete;

    std::string const& name() const { return _name; }
    set_t const& getCids() const { return _cids; }
    void addCids(set_t const& s) { _cids.insert(s.begin(),s.end()); }
    virtual void print( std::ostream& ) const {}
  private:
    set_t _cids;
    const std::string _name;
  };
}

#endif
