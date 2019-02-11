#ifndef DbExample_ConditionEntity2_hh
#define DbExample_ConditionEntity2_hh

#include <set>

namespace mu2e {
  class ConditionsEntity2 {
  public:
    typedef std::shared_ptr<ConditionsEntity2> ptr;
    typedef std::set<int> set_t;

    virtual ~ConditionsEntity2() {}
    virtual std::string const& name() const =0 ;
    set_t const& getCids() const { return _cids; }
    void addCids(set_t const& s) { _cids.insert(s.begin(),s.end()); }
  private:
    set_t _cids;
  };
};

#endif
