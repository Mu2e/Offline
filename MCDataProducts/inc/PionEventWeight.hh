// Sophie Middleton 2021

#ifndef MCDataProducts_inc_PionEventWeight_hh
#define MCDataProducts_inc_PionEventWeight_hh

namespace mu2e {

  class PionEventWeight {
  public:
    double weight() const { return weight_; }

    PionEventWeight(double w) : weight_(w) {}

    // defautl ctr required by ROOT persistency
    PionEventWeight() : weight_(-1) {}

  private:
    double weight_;
  };
}

#endif/*MCDataProducts_inc_PionEventWeight_hh*/
