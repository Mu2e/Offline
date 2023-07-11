#ifndef DbTables_TstCalib2_hh
#define DbTables_TstCalib2_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TstCalib2 : public DbTable {
 public:
  class Row {
   public:
    Row(int channel, std::string status) : _channel(channel), _status(status) {}
    int channel() const { return _channel; }
    std::string const& status() const { return _status; }

   private:
    int _channel;
    std::string _status;
  };

  constexpr static const char* cxname = "TstCalib2";

  TstCalib2() : DbTable(cxname, "tst.calib2", "channel,status") {}

  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  const Row& row(const int channel) const {
    return _rows.at(_chanIndex.at(channel));
  }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  // this is not quite right with strings...
  size_t size() const override {
    size_t b = baseSize() + nrow() * nrow() / 2 + nrow() * sizeof(Row);
    for (auto const& r : _rows) b += r.status().capacity();
    return b;
  };

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoul(columns[0]), columns[1]);
    _chanIndex[_rows.back().channel()] = _rows.size() - 1;
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.channel() << ",";
    sstream << r.status();
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
