#ifndef DbExample_ConditionEntity2_hh
#define DbExample_ConditionEntity2_hh

namespace mu2e {
  class ConditionsEntity2 {
  public:
    typedef std::shared_ptr<const ConditionsEntity2> ptr;
    virtual ~ConditionsEntity2() {}
    virtual std::string const& name() const =0 ;
  };
};

#endif
