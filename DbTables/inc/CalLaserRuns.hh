#ifndef DbTables_CalLaserRuns_hh
#define DbTables_CalLaserRuns_hh

//
// Ad-hoc table (for record-keeping, not calibration itself)
// to record info about a cal laser run
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
    Row(int run, int timeCid, int energyCid, std::string comment) :
      _run(run), _timeCid(timeCid), _energyCid(energyCid), _comment(comment) {}
    int run() const { return _run; }
    int timeCid() const { return _timeCid; }
    int energyCid() const { return _energyCid; }
    std::string comment() const { return _comment; }

   private:
    int _run; // the laser run number
    int _timeCid; // CID of laser run archive table CalLaserTimeCalib
    int _energyCid; // CID of laser run archive table CalLaserEnergyCalib
    std::string _comment;
  };

  constexpr static const char* cxname = "CalLaserRuns";

  CalLaserRuns() : DbTable(cxname, "cal.laserruns", "run,timecid,energycid,comment") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    return baseSize() + nrow() * sizeof(Row);
  };
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stoi(columns[1]),
                       std::stoi(columns[2]), columns[3] );
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.run() << ",";
    sstream << r.timeCid() << ",";
    sstream << r.energyCid() << ",";
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
