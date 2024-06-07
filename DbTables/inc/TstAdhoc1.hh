#ifndef DbTables_TstAdhoc1_hh
#define DbTables_TstAdhoc1_hh

#include "Offline/DbTables/inc/DbTable.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace mu2e {

class TstAdhoc1 : public DbTable {
 public:
  class Row {
   public:
    Row(int exint, float exfloat, std::string exstring) :
        _exint(exint), _exfloat(exfloat), _exstring(exstring) {}
    int exint() const { return _exint; }
    float exfloat() const { return _exfloat; }
    std::string exstring() const { return _exstring; }

   private:
    int _exint;
    float _exfloat;
    std::string _exstring;
  };

  constexpr static const char* cxname = "TstAdhoc1";

  TstAdhoc1() : DbTable(cxname, "tst.adhoc1", "exint,exfloat,exstring") {}
  const Row& rowAt(const std::size_t index) const { return _rows.at(index); }
  std::vector<Row> const& rows() const { return _rows; }
  std::size_t nrow() const override { return _rows.size(); };
  size_t size() const override {
    return baseSize() + nrow() * nrow() / 2 + nrow() * sizeof(Row);
  };
  tableType type() const override { return Adhoc; }

  void addRow(const std::vector<std::string>& columns) override {
    _rows.emplace_back(std::stoi(columns[0]), std::stof(columns[1]),
                       columns[2]);
  }

  void rowToCsv(std::ostringstream& sstream, std::size_t irow) const override {
    Row const& r = _rows.at(irow);
    sstream << r.exint() << ",";
    sstream << std::fixed << std::setprecision(3) << r.exfloat() << ",";
    sstream << r.exstring();
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
