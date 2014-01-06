// Andrei Gaponenko, 2014

#ifndef MCDataProducts_inc_EventWeight_hh
#define MCDataProducts_inc_EventWeight_hh

namespace mu2e {

  class EventWeight {
  public:
    double weight() const { return weight_; }

    EventWeight(double w) : weight_(w) {}

    // defautl ctr required by ROOT persistency
    EventWeight() : weight_(-1) {}

  private:
    double weight_;
  };
}

#endif/*MCDataProducts_inc_EventWeight_hh*/
