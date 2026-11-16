#ifndef DbTables_DAQTiming_hh
#define DbTables_DAQTiming_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <cstdint>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>

namespace mu2e {

class DAQTiming : public DbTable {
 public:
  typedef std::shared_ptr<DAQTiming> ptr_t;
  typedef std::shared_ptr<const DAQTiming> cptr_t;

  class Row {
   public:
    Row( float crvTrackerTimeOffset,
         float caloTrackerTimeOffset,
         float timeFromProtonsToDRMarker) :
        _crvTrackerTimeOffset(crvTrackerTimeOffset),
        _caloTrackerTimeOffset(caloTrackerTimeOffset),
        _timeFromProtonsToDRMarker(timeFromProtonsToDRMarker) {}
    float crvTrackerTimeOffset() const { return _crvTrackerTimeOffset; }
    float caloTrackerTimeOffset() const { return _caloTrackerTimeOffset; }
    float timeFromProtonsToDRMarker() const { return _timeFromProtonsToDRMarker; }

   private:
    float _crvTrackerTimeOffset;
    float _caloTrackerTimeOffset;
    float _timeFromProtonsToDRMarker;
  };

  constexpr static const char* cxname = "DAQTiming";

  DAQTiming() : DbTable(cxname, "daq.timing",
      "crvTrackerTimeOffset,caloTrackerTimeOffset,timeFromProtonsToDRMarker") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  std::size_t nrowFix() const override { return 1; };
  std::size_t size() const override { return baseSize() + nrow() * sizeof(Row); };

  void addRow(const std::vector<std::string>& columns) override {
    if (_rows.size() != 0)
          throw cet::exception("DAQTIMING_BAD_INDEX")
            << "DAQTiming::addRow adding more than one row\n";
    float crvTrackerTimeOffset = std::stof(columns[0]);
    float caloTrackerTimeOffset = std::stof(columns[1]);
    float timeFromProtonsToDRMarker = std::stof(columns[2]);
    _rows.emplace_back(crvTrackerTimeOffset,caloTrackerTimeOffset,timeFromProtonsToDRMarker);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << std::fixed << std::setprecision(3);
    sstream << r.crvTrackerTimeOffset() << ",";
    sstream << r.caloTrackerTimeOffset() << ",";
    sstream << r.timeFromProtonsToDRMarker();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
  }

 private:
  std::vector<Row> _rows;
};

}  // namespace mu2e
#endif
