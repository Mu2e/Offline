#ifndef RecoDataProducts_MSDHit_hh
#define RecoDataProducts_MSDHit_hh

#include <vector>

namespace mu2e {

class MSDHit {
public:
  // constructors
  MSDHit() = default;

  // setters
  void setTime(double time) { _time = time; }
  void setTOT(double tot) { _tot = tot; }

  // accessors
  double time() const { return _time; }
  double tot() const { return _tot; }

  // check if field is present in payload
  bool hasTime() const { return _time >= 0; }
  bool hasTOT() const { return _tot >= 0; }

private:
  // data members
  double _time{-1.0}; // time of hit
  double _tot{-1.0};  // time-over-threshold
};

typedef std::vector<mu2e::MSDHit> MSDHitCollection;

} // namespace mu2e

#endif /* RecoDataProducts_MSDHit_hh */

