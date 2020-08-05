// Andy Edmonds, 2020

#ifndef MCDataProducts_inc_SimStageEfficiency_hh
#define MCDataProducts_inc_SimStageEfficiency_hh

namespace mu2e {

  class SimStageEfficiency {
  public:
    unsigned long numerator() const { return _numerator; }
    unsigned long denominator() const { return _denominator; }

    SimStageEfficiency(unsigned long num, unsigned long denom) : _numerator(num), _denominator(denom) {}

    // defautl ctr required by ROOT persistency
    SimStageEfficiency() : _numerator(0), _denominator(0) {}

    void aggregate(SimStageEfficiency const& other) {
      _numerator += other.numerator();
      _denominator += other.denominator();
    }

    double efficiency() const {
      return (double) _numerator / _denominator;
    }

  private:
    unsigned long _numerator;
    unsigned long _denominator;
  };
}

#endif/*MCDataProducts_inc_SimStageEfficiency_hh*/
