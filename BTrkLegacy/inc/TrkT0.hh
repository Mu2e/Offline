// struct to define t0 as an external parameter
#ifndef TrkT0_hh
#define TrkT0_hh

// simple struct to put together t0 and t0 error
struct TrkT0 {
  TrkT0(double t0, double t0err) : _t0(t0),_t0err(t0err){}
  TrkT0() : _t0(0.0),_t0err(-1.0){}
  double t0() const { return _t0; }
  double t0Err() const { return _t0err; }
  void setT0(double t0, double t0err) { _t0 = t0; _t0err = t0err; }
  double _t0; // estimate of t0 value.  Convention is implementation-specific
  double _t0err; // estimate of t0 error
};

#endif