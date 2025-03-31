#ifndef DbTables_CalLaserRuns_hh
#define DbTables_CalLaserRuns_hh

//
// Ad-hoc table (for record-keeping, not calibration itself)
//

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class CalLaserRuns : public DbTable {
 public:
  class Row {
   public:
    Row(int t0alignRunNumber, int timeCidLaser, std::string comment) :
      _t0alignRunNumber(t0alignRunNumber), _timeCidLaser(timeCidLaser),
      _comment(comment) {}
    int t0alignRunNumber() const { return _t0alignRunNumber; }
    int timeCidLaser() const { return _timeCidLaser; }
    std::string comment() const { return _comment; }

   private:
    int _t0alignRunNumber;
    int _timeCidLaser;
    std::string _comment;
  };

  constexpr static const char* cxname = "CalLaserRuns";

  CalLaserRuns() : DbTable(cxname, "cal.laserruns",
                           "t0alignrunnumber,timecidlaser,comment") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]),
                       columns[2] );
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.t0alignRunNumber() << ",";
    sstream << r.timeCidLaser() << ",";
    sstream << r.comment();
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
