#ifndef DbTables_TstCalib1_hh
#define DbTables_TstCalib1_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TstCalib1 : public DbTable {
 public:
  class Row {
   public:
    Row(int channel, int flag, float dtoe) :
        _channel(channel), _flag(flag), _dtoe(dtoe) {}
    int channel() const { return _channel; }
    int flag() const { return _flag; }
    float dToE() const { return _dtoe; }

   private:
    int _channel;
    int _flag;
    float _dtoe;
  };

  constexpr static const char* cxname = "TstCalib1";

  TstCalib1() : DbTable(cxname, "tst.calib1", "channel,flag,dtoe") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(const int channel) const {
    return _rows.at(_chanIndex.at(channel));
  }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  // this table should always be 3 rows
  virtual std::size_t nrowFix() const override { return 3; };
  size_t size() const override {
    return baseSize() + nrow() * nrow() / 2 + nrow() * sizeof(Row);
  };

  void addRow(const std::vector<std::string>& columns) override {
    int channel = std::stoi(columns[0]);
    // enforce a strict sequential order - optional
    if (channel != int(_rows.size())) {
      throw cet::exception("TSTCALIB1_BAD_INDEX")
          << "TstCalib1::addRow found index out of order: " << channel
          << " != " << _rows.back().channel() + 1 << "\n";
    }
    _rows.emplace_back(channel, std::stoi(columns[1]), std::stof(columns[2]));
    // add this channel to the map index - optional
    _chanIndex[_rows.back().channel()] = _rows.size() - 1;
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << r.flag() << ",";
    sstream << std::fixed << std::setprecision(3) << r.dToE();
  }

  virtual void clear() override {
    baseClear();
    _rows.clear();
    _chanIndex.clear();
  }

 private:
  std::vector<Row> _rows;
  std::map<int, std::size_t> _chanIndex;
};

}  // namespace mu2e
#endif
