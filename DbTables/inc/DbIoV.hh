#ifndef DbTables_DbIoV_hh
#define DbTables_DbIoV_hh

#include <climits>
#include <cstdint>
#include <iomanip>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace mu2e {

class DbIoV {
 public:
  DbIoV() { setNull(); }  // defaults to no valid range
  DbIoV(uint32_t startRun, uint32_t startSubrun, uint32_t endRun,
        uint32_t endSubrun) :
      _startRun(startRun),
      _startSubrun(startSubrun), _endRun(endRun), _endSubrun(endSubrun) {}
  DbIoV(std::string const& text) {
    setByString(text);
  }  // defaults to no valid range

  void set(uint32_t startRun, uint32_t startSubrun, uint32_t endRun,
           uint32_t endSubrun) {
    _startRun = startRun;
    _startSubrun = startSubrun;
    _endRun = endRun;
    _endSubrun = endSubrun;
  }

  void setMax() {
    _startRun = 0;
    _startSubrun = 0;
    _endRun = maxRun();
    _endSubrun = maxSubrun();
  }

  void setNull() {
    _startRun = 0;
    _startSubrun = 0;
    _endRun = 0;
    _endSubrun = 0;
  }

  // remove argument's interval from this
  void subtract(DbIoV const& iov, uint32_t run = 0, uint32_t subrun = 0);
  // remove any interval not overlapping
  void overlap(DbIoV const& iov);

  // formats (see the wiki for more info):
  // start_run:start_subrun-end_run:end_subrun
  // start_run:start_subrun-end_run
  // start_run-end_run
  // start_run
  void setByString(std::string s);

  bool inInterval(uint32_t run, uint32_t subrun) const {
    if (run < _startRun) return false;
    if (run == _startRun && subrun < _startSubrun) return false;
    if (run > _endRun) return false;
    if (run == _endRun && subrun > _endSubrun) return false;
    return true;
  }

  // 0 not overlapping, 1 complete, 2 non-overlapping includes start only,
  // 3 non-overlapping include end only, 4 splits the iov
  int isOverlapping(DbIoV const& iov) const;

  bool isNull() const {
    if (_startRun == 0 && _startSubrun == 0 && _endRun == 0 && _endSubrun == 0)
      return true;
    if (_startRun > _endRun) return true;
    if (_startRun == _endRun && _startSubrun > _endSubrun) return true;
    return false;
  }

  uint32_t startRun() const { return _startRun; }
  uint32_t startSubrun() const { return _startSubrun; }
  uint32_t endRun() const { return _endRun; }
  uint32_t endSubrun() const { return _endSubrun; }

  constexpr static uint32_t maxRun() { return 999999; }
  constexpr static uint32_t maxSubrun() { return 999999; }

  // in format: run subrun run subrun
  std::string simpleString() const;
  // in format run:subrun-run:subrun, with fixed width or compressed
  std::string to_string(bool compress = false) const;

 private:
  uint32_t _startRun;
  uint32_t _startSubrun;
  uint32_t _endRun;
  uint32_t _endSubrun;
};

inline std::ostream& operator<<(std::ostream& ost, mu2e::DbIoV const& iov) {
  ost << iov.to_string(true);
  return ost;
}

}  // namespace mu2e
#endif
