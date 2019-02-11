#ifndef DbTables_DbIoV_hh
#define DbTables_DbIoV_hh

#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <ostream>
#include <iomanip>
#include <cstdint>
#include <climits>

namespace mu2e {

  class DbIoV {
  public:

    DbIoV() { set(0,0,0,0); } // defaults to no valid range
    DbIoV(uint32_t startRun, uint32_t startSubrun,
	  uint32_t endRun, uint32_t endSubrun):
      _startRun(startRun),_startSubrun(startSubrun),
      _endRun(endRun),_endSubrun(endSubrun) {} 

    void set(uint32_t startRun, uint32_t startSubrun,
	     uint32_t endRun, uint32_t endSubrun) {
      _startRun = startRun;
      _startSubrun = startSubrun;
      _endRun = endRun;
      _endSubrun = endSubrun;
    }

    void setMax() {
      _startRun = 0;
      _startSubrun = 0;
      _endRun = maxr;
      _endSubrun = maxsr;
    }

    // remove argument's interval from this
    void subtract(DbIoV const& iov, uint32_t run=0, uint32_t subrun=0);
    // remove any interval not overlapping
    void overlap(DbIoV const& iov);

    // formats (see the wiki for more info):
    // start_run:start_subrun-end_run:end_subrun
    // start_run:start_subrun-end_subrun
    // start_run-end_subrun
    // start_run
    void setByString(std::string s);

    bool inInterval(uint32_t run, uint32_t subrun) const {
      if(run<_startRun) return false;
      if(run==_startRun and subrun<_startSubrun) return false;
      if(run>_endRun) return false;
      if(run==_endRun and subrun>_endSubrun) return false;
      return true;
    }

    uint32_t startRun() const {return _startRun;}
    uint32_t startSubrun() const {return _startSubrun;}
    uint32_t endRun() const {return _endRun;}
    uint32_t endSubrun() const {return _endSubrun;}

    uint32_t maxRun() { return maxr; }
    uint32_t maxSubrun() { return maxsr; }

    // in format: run subrun run subrun
    std::string simpleString() const;
    // in format run:subrun-run:subrun, with fixed width or compressed
    std::string to_string(bool compress = false) const;

  private:
    static const uint32_t maxsr = 999999;
    static const uint32_t maxr  = 999999;
    
    uint32_t _startRun;
    uint32_t _startSubrun;
    uint32_t _endRun;
    uint32_t _endSubrun;
  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   mu2e::DbIoV const& iov) {
    ost << iov.to_string(true);
    return ost;
  }

}
#endif
