// Michael MacKenzie, 2025

#ifndef MCDataProducts_inc_SumOfWeights_hh
#define MCDataProducts_inc_SumOfWeights_hh

namespace mu2e {

  class SumOfWeights {
  public:

    // constructors
    SumOfWeights() : sum_(0.), count_(0) {}
    SumOfWeights(double sum, unsigned long count) : sum_(sum), count_(count) {}

    // accessors
    double sum  () const { return sum_  ; }
    unsigned long   count() const { return count_; }

    // functions
    void add(const double weight, const unsigned long count = 1) {
      sum_ += weight;
      count_ += count;
    }
    void reset() {
      sum_ = 0.;
      count_ = 0;
    }
    double avg() const {
      return (count_ > 0) ? sum_ / count_ : 0.;
    }

  private:
    double sum_  ; // total sum of weights
    unsigned long   count_; // N(events) in the sum
  };
}

#endif/*MCDataProducts_inc_EventWeight_hh*/
