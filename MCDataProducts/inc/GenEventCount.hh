// An object of GenEventCount type is used to record the number of
// events generated in a job, before any filters are applied.
// Production is structured to have a separate subrun for each first
// stage job, and GenEventCount is recorded in the SubRun object.
//
// Andrei Gaponenko, 2013

#ifndef MCDataProducts_inc_GenEventCount_hh
#define MCDataProducts_inc_GenEventCount_hh

namespace mu2e {

  class GenEventCount {
  public:
    typedef unsigned long count_t;

    count_t count() const { return count_; }

    GenEventCount(unsigned long n) : count_(n) {}

    // defautl ctr required by ROOT persistency
    GenEventCount() : count_(-1) {}

  private:
    count_t count_;
  };
}

#endif/*MCDataProducts_inc_GenEventCount_hh*/
